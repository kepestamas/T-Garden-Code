import pandas as pd
import graph
import math
from itertools import combinations
import networkx as nx
import random
import os
import glob
import multiprocessing
import geology_costs

import numpy as np


from shapely.geometry import LineString
from pyproj import Geod

geod = Geod(ellps="WGS84")

def geodesic_line(lon1, lat1, lon2, lat2, n_points=200):
    all_points = geod.npts(lon1, lat1, lon2, lat2, npts=n_points)
    lons, lats = zip(*all_points)
    lons = list(lons)
    lats = list(lats)
    return LineString(zip(lons, lats))

lines = {}

def arcs_intersect(p1, p2, q1, q2):
    if str(p1) not in lines.keys():
        lines[str(p1)] = {}
    if str(p2) not in lines[str(p1)].keys():
        lines[str(p1)][str(p2)] = geodesic_line(p1[1], p1[0], p2[1], p2[0]) #p and q are given az two points with (lat, lon) encoding
    line1 = lines[str(p1)][str(p2)]

    if str(q1) not in lines.keys():
        lines[str(q1)] = {}
    if str(q2) not in lines[str(q1)].keys():
        lines[str(q1)][str(q2)] = geodesic_line(q1[1], q1[0], q2[1], q2[0]) #p and q are given az two points with (lat, lon) encoding
    line2 = lines[str(q1)][str(q2)]


    return line1.intersects(line2)




def populate_cities_dictionary():
    df = pd.read_csv('../inputs/geonames/cities1000.txt', sep='\t', header=None, names=[
        'geonameid', 'name', 'asciiname', 'alternatenames', 'latitude', 'longitude',
        'feature class', 'feature code', 'country code', 'cc2', 'admin1 code',
        'admin2 code', 'admin3 code', 'admin4 code', 'population', 'elevation',
        'dem', 'timezone', 'modification date'
    ])

    regions = {}
    country_df = pd.read_csv('../inputs/geonames/countries.txt', sep='\t', header=None)
    countries_by_code = {}
    for index, code  in enumerate(country_df[0]):
        countries_by_code[code] = country_df[4][index]
    cities = []
    for index, city_name  in enumerate(df["asciiname"]):
        city = {}
        city["name"] = city_name
        city["longitude"] = df["longitude"][index]
        city["latitude"] = df["latitude"][index]
        city["region"] = df["timezone"][index].split("/")[0]
        city["population"] = df["population"][index]
        city["country"] = countries_by_code[df['country code'][index]]
        city["country code"] = df['country code'][index]
        cities.append(city)
    cities_by_regions = {}
    regions = []


    for city in cities:
        if city["country code"] in ["US"]:
            if "US" not in cities_by_regions.keys():
                regions.append("US")
                cities_by_regions["US"] = [city]
            else:
                cities_by_regions["US"].append(city)
        if city["country code"] in ["US", "CA"]:
            city["region"] = "North America"
        elif city["country code"] in ["BR", "AR", "CL", "BO", "PE", "UY", "PY", "EC", "CO",  "VE", "GY", "SR", "GF", "MX"]:
            city["region"] = "South America"
        elif city["region"] == "America":
            city["region"] = "South America"
        elif city["country code"] in ["MA", "LY", "TN", "DZ", "EG", "EH", "SD"]:
            city["region"] = "North Africa"
        elif city["region"] == "Africa":
            city["region"] = "Sub-Saharan Africa"
        elif city["country code"] in ["SY", "LB", "IQ", "IR", "IL", "JO", "PS", "OM", "YE", "SA", "AE", "BH", "QA", "KW", "AZ", "AM", "GE", "CY"]:
            city["region"] = "Middle East"
        elif city["country code"] in ["JP", "TW", "CN", "KR", "KP", "MN", "MO", "HK"]:
            city["region"] = "East Asia"
        elif city["country code"] in ["TJ", "TM", "UZ", "AF", "KG", "KZ"]:
            city["region"] = "Central Asia"
        elif city["country code"] in ["PK", "NP", "IN", "BD", "BT", "MV", "LK", "MM" ]:
            city["region"] = "South Asia"
        elif city["country code"] in ["RU"] and city["region"] == "Asia":
            city["region"] = "Siberia"
        elif city["region"] == "Asia":
            city["region"] = "South-East Asia"  
        if city["region"] not in cities_by_regions.keys():
            regions.append(city["region"])
            cities_by_regions[city["region"]] = [city]
        else:
            cities_by_regions[city["region"]].append(city)

    for region in regions:
        cities_by_regions[region] = sorted(cities_by_regions[region], key=lambda d: d["population"], reverse=True)

    return cities_by_regions


#  Haversine formula for distance
def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0  # Earth radius in kilometers
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2)**2
    return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

#TODO cost calculation here
def calculate_costs(distances, cities_dict):
    print("Calculating costs")
    costs = distances.copy()
    for city_pair in costs.keys():
        # print(city_pair)
        city_a = cities_dict[city_pair[0]]
        city_b = cities_dict[city_pair[1]]
        if city_a["country"] == city_b["country"]:
            # print("Same country, lowering cost " + str(city_a["name"]) + " " + str(city_b["name"]))
            costs[city_pair] *= 0.5 # Half cost for same country connections
        if geology_costs.crosses_water(city_a["longitude"], city_a["latitude"], city_b["longitude"], city_b["latitude"], 30000):
            # print("Crosses water, lowering cost " + str(city_a["name"]) + " " + str(city_b["name"]))
            costs[city_pair] *= 0.7 # 70% cost for water crossing
    return costs

# Generate graph with realistic edges
def generate_edges(cities, edge_count, initial_edges = None, shuffle_edges = False):
    distances = {
        (a['name'], b['name']): haversine(a["latitude"], a["longitude"], b["latitude"], b["longitude"])
        for a, b in combinations(cities, 2)
    }
    cities_dict = {c["name"]: c for c in cities}
    
    
    if initial_edges == None:
        edges = generate_complete_graph_edges(edge_count, distances)
    else:
        costs = calculate_costs(distances, cities_dict) 
        edge_count -= len(initial_edges) # the mst edges should also count toward edge count
        
        sorted_edges = sorted(costs.items(), key=lambda x: x[1])
        edges = convert_initial_edges_to_usable_edges(cities, initial_edges)
        sorted_edges = filter_based_on_edges(edges, cities_dict, sorted_edges)
        print("Added initial edges, filtered all remaining edges to remove intersects")
        edge_i = 0
        i = 0
        global REFUSAL_RATE
        while True:
            (u, v), dist = sorted_edges[0]
            should_refuse = random.random()
            while should_refuse < REFUSAL_RATE and len(sorted_edges) > edge_count - edge_i: # if there aren't enough edges, then no longer consider refusal
                sorted_edges = sorted_edges[1:]
                (u, v), dist = sorted_edges[0]
                should_refuse = random.random() + (REFUSAL_RATE / 2) #! added a +0.1 so that the second edge, has a larger chance of being selected, so that this will not fail
            print(f'Considering edges no. {i} out of {len(sorted_edges)}')
            # if (u, v) not in initial_edges and (v, u) not in initial_edges: #! Placed the filter at the distance calculator
            was = False
            if not was:
                edges.append([u, v, dist])
                sorted_edges = filter_based_on_edges([[u, v, dist]], cities_dict, sorted_edges)
                edge_i += 1
                if edge_i % 1 == 0:
                    print(f'Found edge no. {edge_i}, remaining {edge_count - edge_i}')
                if edge_i == edge_count:
                    print(f'Final edge found')
                    break
            i += 1
            if len(sorted_edges) == 0:
                print("Finished considering all edges")
                break
    return edges

def filter_based_on_edges(edges, cities_dict, sorted_edges):
    i = 0
    for edge in edges:
        print(f'Filtering based on edges no. {i} out of the given edges, {len(sorted_edges)} edges remain')
        i += 1
        new_sorted = []
        for sorted_edge in sorted_edges:
            if not arcs_intersect(
                (cities_dict[edge[0]]["latitude"], cities_dict[edge[0]]["longitude"]),
                (cities_dict[edge[1]]["latitude"], cities_dict[edge[1]]["longitude"]), 
                (cities_dict[sorted_edge[0][0]]["latitude"], cities_dict[sorted_edge[0][0]]["longitude"]), 
                (cities_dict[sorted_edge[0][1]]["latitude"], cities_dict[sorted_edge[0][1]]["longitude"])):
                new_sorted.append(sorted_edge)
        sorted_edges = new_sorted.copy()
            
    return sorted_edges

def convert_initial_edges_to_usable_edges(cities, initial_edges):
    edges = []
    for edge in initial_edges:
        for city in cities:
            if city['name'] == edge[0]:
                a = city
            elif city['name'] == edge[1]:
                b = city
        initial_edge_with_distance = [a['name'], b['name'], haversine(a["latitude"], a["longitude"], b["latitude"], b["longitude"])]
        edges.append(initial_edge_with_distance)
    return edges

#TODO: Look at maybe not needing the check, since we need all sorted edges.
def generate_complete_graph_edges(edge_count, distances):
    print("Generating all edges")
    sorted_edges = sorted(distances.items(), key=lambda x: x[1])
    edges = []
    edge_i = 0
    for (u, v), dist in sorted_edges:
        # if (u, v) not in edges and (v, u) not in edges:
        edges.append([u, v, dist])
        edge_i += 1
        if edge_i == edge_count:
            break
    return edges


def print_network_to_file(path, cities, edges):
    f = open(path, "w+")
    for city in cities:
        city_str = str(city["name"]) + "|" + str(city["longitude"]) + "|" + str(city["latitude"])
        print(city_str, file=f)
    for edge in edges:
        if len(edge) > 2:
            dist = edge[2]
        else:
            dist = 0
        edge_str = str(edge[0]) + "|" + str(edge[1]) + "|" + str(dist)
        print(edge_str, file=f)
    f.close()


def generate_network(cities_by_regions: dict, region = "Europe", node_count = 100, edge_count = 200, randomize_selection = False, mst=2):
    global PMST_RATE
    file_name = generate_new_network_file_name(region, node_count, edge_count, randomize_selection, mst)
    folder_name = "../networks/generated/"
    path =  folder_name + file_name
    tmp_path = path + ".tmp"
    possible_cities = select_possible_cities_from_region_data(cities_by_regions, region)
    selected_cities = select_cities_for_network(node_count, randomize_selection, possible_cities)
    if len(selected_cities) != node_count:
        print("!!!!!!!!!!!!!!!City selection failed!!!!!!!!!!!!!!!!!!!!!!")
    if mst > 0:
        full_edges = generate_edges(selected_cities, node_count * node_count)
        initial_edges = extract_initial_edges_using_mst(file_name, tmp_path, selected_cities, full_edges)
        if mst == 1:
            initial_edges = random.sample(initial_edges, int(len(initial_edges) * PMST_RATE))
    else:
        initial_edges = []
    final_edges = generate_edges(selected_cities, edge_count=edge_count, initial_edges=initial_edges)
    if mst < 2: 
        final_edges = ensure_connection(file_name, tmp_path, selected_cities, final_edges)
    print_network_to_file(path, selected_cities, final_edges)


def ensure_connection(file_name, tmp_path, selected_cities, initial_final_edges):
    print_network_to_file(tmp_path, selected_cities, initial_final_edges)
    print("Ensuring connected graph")
    G = graph.Graph(source_name="generated", file_name =  file_name + ".tmp", format="generated")
    # G.G = nx.minimum_spanning_tree(G.G) #! Ensure connection here
    
    nr_added = 0
    while not nx.is_connected(G.G):
        components = list(nx.connected_components(G.G))
        best_dist = float("inf")
        best_edge = None
        for i, comp1 in enumerate(components):
            for j, comp2 in enumerate(components):
                if i >= j:
                    continue
                for u in comp1:
                    lon1, lat1 = G.G.nodes[u]["Longitude"], G.G.nodes[u]["Latitude"]
                    for v in comp2:
                        lon2, lat2 = G.G.nodes[v]["Longitude"], G.G.nodes[v]["Latitude"]
                        d = haversine(lat1, lon1, lat2, lon2)
                        if d < best_dist:
                            best_dist = d
                            best_edge = (u, v)
        u, v = best_edge
        G.G.add_edge(u, v)
        nr_added += 1

    edges = list(G.G.edges())
    random.shuffle(edges)
    nr_removed = 0
    for u, v in edges:
        if nr_removed >= nr_added:
            break
        if (u, v) not in nx.bridges(G.G):  # only remove if not a bridge
            G.G.remove_edge(u, v)
            nr_removed += 1

    folder = "../networks/generated"
    extension = ".tmp"  
    for file_path in glob.glob(os.path.join(folder, f"*{extension}")):
        os.remove(file_path)
        print(f"Deleted: {file_path}")

    final_edges_with_dist = []
    for final_edge in G.G.edges():
        a = G.G.nodes[final_edge[0]]
        b = G.G.nodes[final_edge[1]]
        final_edges_with_dist.append((final_edge[0], final_edge[1], haversine(a["Latitude"], a["Longitude"], b["Latitude"], b["Longitude"])))
    return final_edges_with_dist


def extract_initial_edges_using_mst(file_name, tmp_path, selected_cities, full_edges):
    initial_edges = []
    print_network_to_file(tmp_path, selected_cities, full_edges)
    print("Generating mst")
    G = graph.Graph(source_name="generated", file_name =  file_name + ".tmp", format="generated")
    G.G = nx.minimum_spanning_tree(G.G)
    for edge in G.G.edges():
        initial_edges.append(edge)
    folder = "../networks/generated"
    extension = ".tmp"  
    for file_path in glob.glob(os.path.join(folder, f"*{extension}")):
        os.remove(file_path)
        print(f"Deleted: {file_path}")
    return initial_edges


def too_close(city, selected_cities):
    for selected_city in selected_cities:
        _, _, dist = geod.inv(city["longitude"], city["latitude"], selected_city["longitude"], selected_city["latitude"])
        if dist/1000 < CLOSENESS_LIMIT:
            return True
    return False


def select_cities_for_network(node_count, randomize_selection, possible_cities):
    selected_cities = []
    for i in range(node_count):
        if not randomize_selection:
            city = possible_cities[i] 
        else:
            weights = [1 / (i + 1) for i in range(len(possible_cities))]
            city = random.choices(possible_cities, weights=weights, k = 1)[0]
            tries = 0
            while too_close(city, selected_cities) and tries < node_count:
                city = random.choices(possible_cities, weights=weights, k = 1)[0]
                tries += 1
        selected_cities.append(city)
    return selected_cities

def select_possible_cities_from_region_data(cities_by_regions, region):
    possible_cities = []
    if region != "Global":
        regions = region.split("+")
        for region_it in regions:
            if region_it == "America":
                possible_cities += cities_by_regions["North America"]
                # possible_cities += cities_by_regions["Central America and Carribean"]
                possible_cities += cities_by_regions["South America"]
            elif region_it == "Africa":
                possible_cities += cities_by_regions["Sub-Saharan Africa"]
                possible_cities += cities_by_regions["North Africa"]
            elif region_it == "Middle East":
                possible_cities += cities_by_regions["Middle East"]
                for city in cities_by_regions["Europe"]:
                    if city["country code"] in ["TR"]:
                        possible_cities.append(city)
                for city in cities_by_regions["North Africa"]:
                    if city["country code"] in ["EG"]:
                        possible_cities.append(city)       
            elif region_it == "Asia":
                possible_cities += cities_by_regions["Middle East"]
                possible_cities += cities_by_regions["Central Asia"]
                possible_cities += cities_by_regions["Siberia"]
                possible_cities += cities_by_regions["South Asia"]
                possible_cities += cities_by_regions["South-East Asia"]
                possible_cities += cities_by_regions["East Asia"]
            else:
                possible_cities += cities_by_regions[region_it]
    else:
        for region_it in cities_by_regions.keys():
            possible_cities += cities_by_regions[region_it]
    possible_cities = sorted(possible_cities, key=lambda d: d["population"], reverse=True)
    return possible_cities

def generate_new_network_file_name(region, node_count, edge_count, randomize_selection, mst):
    file_name = region  + "_" + str(node_count) + "_" + str(edge_count)
    
    mst_str = "nmst"
    if mst == 1:
        mst_str = "pmst"
    elif mst == 2:
        mst_str = "mst"
    file_name +=  "_" + mst_str

    if randomize_selection:
        file_name += "_rand"
    file_name += ".txt"

    return file_name



# region = "Asia+Europe+Africa"
# node_count = 100
# randomize_selection = True
# mst = 1 # 0 - no mst, 1 - partial mst 50%, 2 - full mst

CLOSENESS_LIMIT = 50 # kms, only applies to random
REFUSAL_RATE = 0.1
PMST_RATE = 0.75


def generate_networks_from_parameters(regions, node_edge_counts, randomize_selections, msts):
    # cities_by_regions = populate_cities_dictionary()

    for region in regions:
        for node_count, edge_count in node_edge_counts:
            for randomize_selection in randomize_selections:
                for mst in msts:
                    # generate_network(cities_by_regions, region=region, node_count=node_count, edge_count=edge_count, randomize_selection=randomize_selection, mst=mst)

                    random_tag = ""
                    if randomize_selection:
                        random_tag = "_rand"

                    mst_str = "_nmst"
                    if mst == 1:
                        mst_str = "_pmst"
                    elif mst == 2:
                        mst_str = "_mst"

                    # G = graph.Graph(source_name="generated", file_name =  region + "_" + str(node_count) + "_" + str(edge_count) + mst_str +  random_tag + ".txt", format="generated")
                    G = graph.Graph(source_name="unified", file_name =  region + "_" + str(node_count) + "_" + str(edge_count) + mst_str +  random_tag + ".gml", format="unified")
                    # print(f"Looking at {G.G.graph["Network"]}")
                    # G.remove_triangular_edges()
                    # G.serialize_to_gml()
                    # G.show_map_plot()
                    G.save_map_plot()




