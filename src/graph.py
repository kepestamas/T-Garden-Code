import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import requests
from pyproj import Geod


class Graph:
    def __init__(self, source_name, file_name, format = "gml", ignore_atlantic = True):

        self.modified_by_drawing = False
        self.drawing_selected_nodes = []
        self.drawing_selected_edge = None
        self.source_name = source_name
        self.file_name = file_name
        self.format = format
        self.path = "../networks/" + source_name + "/" + file_name
        if source_name == "unified":
            self.G : nx.Graph = nx.read_gml(self.path, label="id")
            self.source_name = self.G.graph["Source"]
        else:
            if format == "gml":
                self.G : nx.Graph = nx.read_gml(self.path, label="id")
            if format == "sndlib":
                self.G : nx.Graph = self.__read_sndlib(self.path)
            if format == "generated":
                self.G : nx.Graph = self.__read_generated(self.path, file_name)
        #! Removed for generator testing
        if source_name == "topology_zoo":
            self.__remove_non_internal_nodes() 
            self.__remove_non_connected_nodes()
            self.__get_city_coords_for_non_descriptive_nodes()
            self.__calculate_hyperedge_locations()
            self.__calculate_hyperedge_locations() # a second calculation for more accurate hyperedges
        
        self.should_route = False
        if not ignore_atlantic:
            self._prepare_atlantic_traversal_if_necessary()
        self.latitude_limits = (-360, 360)
        self.longitude_limits = (-360, 360)
        self.geod = Geod(ellps="WGS84")

    def __read_generated(self, path, file_name):
        G : nx.Graph = nx.Graph()
        f = open(path, "r")
        file_name_params = file_name.split("_")
        G.graph["Network"] = file_name.replace(".txt", "")
        G.graph["GeoLocation"] = file_name_params[0]
        G.graph["GeoExtent"] = "Continent+"
        node_count = int(file_name_params[1].split(".")[0])
        edge_i = 0
        lines = f.readlines()
        for index, line in enumerate(lines):
            split_line = line.strip().split("|")
            if index < node_count:
                G.add_node(split_line[0])
                nx.set_node_attributes(G, {split_line[0]: {"Longitude" : float(split_line[1]), "Latitude" : float(split_line[2]), "label" : split_line[0], "Internal" : 1}})
            else:
                G.add_edge(split_line[0], split_line[1])
                nx.set_edge_attributes(G, {(split_line[0], split_line[1]) : {"id": "E" + str(edge_i), "weight": float(split_line[2])}})
                edge_i += 1
        return G


    def serialize_to_gml(self):
        print("Serializing")
        self.__remove_non_connected_nodes() # added this, for extremal node removal in editing process
        modified_str = ""
        if self.modified_by_drawing:
            modified_str = "modified/"
        f = open(("../networks/unified/" + modified_str + self.file_name).replace(".txt",".gml"), "w+")
        content = "graph [\n" + self.__meta_to_gml() + self.__nodes_to_gml() + self.__edges_to_gml() + "]"
        print(content, file=f)
        f.close()

    def edge_length_in_km(self, edge):
        _, _, dist = self.geod.inv(self.G.nodes[edge[0]]["Longitude"], self.G.nodes[edge[0]]["Latitude"], self.G.nodes[edge[1]]["Longitude"], self.G.nodes[edge[1]]["Latitude"])
        return dist / 1000

    def remove_triangular_edges(self):
        print(f"Edge count: {len(list(self.G.edges))}")
        for node1_id in self.G.nodes:
            for node2_id in self.G.nodes:
                for node3_id in self.G.nodes:
                    if (node1_id, node2_id) in self.G.edges and (node2_id, node3_id) in self.G.edges and (node1_id, node3_id) in self.G.edges:
                        _, _, distance1 = self.geod.inv(self.G.nodes[node1_id]["Longitude"], self.G.nodes[node1_id]["Latitude"], self.G.nodes[node2_id]["Longitude"], self.G.nodes[node2_id]["Latitude"])
                        _, _, distance2 = self.geod.inv(self.G.nodes[node2_id]["Longitude"], self.G.nodes[node2_id]["Latitude"], self.G.nodes[node3_id]["Longitude"], self.G.nodes[node3_id]["Latitude"])
                        _, _, distance3 = self.geod.inv(self.G.nodes[node3_id]["Longitude"], self.G.nodes[node3_id]["Latitude"], self.G.nodes[node1_id]["Longitude"], self.G.nodes[node1_id]["Latitude"])
                        old_distance = [distance1, distance2, distance3]
                        if distance1 > distance2:
                            distance1, distance2 = distance2, distance1
                        if distance2 > distance3:
                            distance3, distance2 = distance2, distance3
                        if distance1 > distance2:
                            distance1, distance2 = distance2, distance1
                        removed_edge = None
                        if distance3 >= 0.9 * (distance1 + distance2):
                            if distance3 == old_distance[0]:
                                self.G.remove_edge(node1_id, node2_id)
                                removed_edge = (node1_id, node2_id)
                            elif distance3 == old_distance[1]:
                                self.G.remove_edge(node2_id, node3_id)
                                removed_edge = (node2_id, node3_id)
                            else:
                                self.G.remove_edge(node1_id, node3_id)
                                removed_edge = (node1_id, node3_id)
                            if not nx.is_connected(self.G):
                                self.G.add_edge(removed_edge[0], removed_edge[1])
        print(f"Final Edge count: {len(list(self.G.edges))}")

    def __meta_to_gml(self):
        prefix = "  "
        data = ""
        data += prefix + "multigraph 1\n"
        data += prefix + "Network \"" + self.G.graph["Network"].strip() + "\"\n"
        data += prefix + "Source \"" + self.source_name + "\"\n"
        data += prefix + "SourceFormat \"" + self.format + "\"\n"
        data += prefix + "GeoExtent \"" + self.G.graph["GeoExtent"].strip() + "\"\n"
        data += prefix + "GeoLocation \"" + self.G.graph["GeoLocation"].strip() + "\"\n"
        return data
    
    def __nodes_to_gml(self):
        prefix = "  "
        data = ""

        hyperedge_count = 0
        for node_id in self.G.nodes:
            data += prefix + "node [\n"
            # print(self.G.nodes[node_id])
            data += prefix + prefix + "id \"" + str(node_id) + "\"\n"
            if "hyperedge" not in self.G.nodes[node_id].keys():
                data += prefix + prefix + "label \"" + str(self.G.nodes[node_id]["label"]) + "\"\n"
            else:
                data += prefix + prefix + "label \"" + "Hyperedge_" + str(hyperedge_count) + "\"\n"
                hyperedge_count += 1
                data += prefix + prefix + "hyperedge 1\n"
            data += prefix + prefix + "Internal 1\n"
            data += prefix + prefix + "Longitude " + str(self.G.nodes[node_id]["Longitude"]) + "\n"
            data += prefix + prefix + "Latitude " + str(self.G.nodes[node_id]["Latitude"]) + "\n"    
            data += prefix + "]\n"
        return data
    
    def __edges_to_gml(self):
        prefix = "  "
        data = ""
        labelless_count = 0
        for edge in self.G.edges:
            data += prefix + "edge [\n"
            data += prefix + prefix + "source \"" + str(edge[0]) + "\"\n"
            data += prefix + prefix + "target \"" + str(edge[1]) + "\"\n"
            if "id" in self.G.edges[edge].keys():
                data += prefix + prefix + "id \"" + str(self.G.edges[edge]["id"]) + "\"\n"
            else:
                data += prefix + prefix + "id \"" + "Non_labeled_" + str(labelless_count) + "\"\n"
                labelless_count += 1
            data += "]\n"
        return data

    def __get_city_coords_for_non_descriptive_nodes(self):
        for node_id in self.G.nodes:
            if "Latitude" not in self.G.nodes[node_id].keys() and "hyperedge" not in self.G.nodes[node_id].keys() and ("label" in self.G.nodes[node_id].keys() and self.G.nodes[node_id]["label"] != "None"):
                address = str(self.G.nodes[node_id]["label"])
                print(address)
                headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_5) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/50.0.2661.102 Safari/537.36'}
                url = 'https://nominatim.openstreetmap.org/search?q=' + address +'&format=json'

                response = requests.get(url, headers=headers).json()
                # print(response)
                
                nx.set_node_attributes(self.G, {node_id: {"Longitude" : float(response[0]["lon"]), "Latitude" : float(response[0]["lat"])}})

    def __remove_non_connected_nodes(self):
        non_connected_nodes = []
        for node_id in self.G.nodes:
            non_connected = True
            for edge in self.G.edges:
                if node_id == edge[0] or node_id == edge[1]:
                    non_connected = False
                    break
            if non_connected:
                non_connected_nodes.append(node_id)
        self.G.remove_nodes_from(non_connected_nodes)

    def __remove_non_internal_nodes(self):
        non_internal_nodes = []
        for node_id in self.G.nodes:
            if "Internal" in self.G.nodes[node_id].keys() and self.G.nodes[node_id]["Internal"] != 1:
                non_internal_nodes.append(node_id)
        self.G.remove_nodes_from(non_internal_nodes)

    def __read_sndlib(self, path):
        G : nx.Graph = nx.Graph()
        f = open(path, "r")
        current_section = ""
        for line in f.readlines():
            line = line.strip()
            if line == ")":
                current_section = ""
            if "=" in line and current_section == "meta":
                split_line = line.split("=")
                if "name" in split_line[0]:
                    G.graph["Network"] = split_line[1].strip()
                if "geoextend" in split_line[0]:
                    G.graph["GeoExtent"] = split_line[1].strip()
                if "geolocation" in split_line[0]:
                    G.graph["GeoLocation"] = split_line[1].strip()
            if current_section == "nodes":
                split_line = line.split(" ")
                G.add_node(split_line[0])
                nx.set_node_attributes(G, {split_line[0]: {"Longitude" : float(split_line[2]), "Latitude" : float(split_line[3]), "label" : split_line[0], "Internal" : 1}})
            if current_section == "edges":
                split_line = line.split(" ")
                G.add_edge(split_line[2], split_line[3])
                nx.set_edge_attributes(G, {(split_line[2], split_line[3]) : {"id": split_line[0]}})
            if "META (" in line:
                current_section = "meta"
            if "NODES (" in line:
                current_section = "nodes"
            if "LINKS (" in line:
                current_section = "edges"
        f.close()
        return G
    
    def get_node(self, node_id):
        return self.G.nodes[node_id]

    def get_nodes(self):
        return self.G.nodes
    
    def get_network_info_dictionary(self):
        network_dict = {}
        network_dict["source"] = self.source_name
        network_dict["name"] = self.G.graph["Network"]
        network_dict["geo_extent"] = self.G.graph["GeoExtent"]
        network_dict["geo_location"] =  self.G.graph["GeoLocation"]
        network_dict["node_count"] = len(self.G.nodes)
        network_dict["edge_count"] = len(self.G.edges)
        network_dict["download_link"] = self.path.replace("../", "/")
        return network_dict

    def get_latitude_list(self):
        latitudes = []
        for node_id in self.G.nodes:
            if "Latitude" in self.G.nodes[node_id].keys():
                latitudes.append(self.G.nodes[node_id]["Latitude"])
        return latitudes

    def get_hyperedge_latitude_list(self):
        latitudes = []
        for node_id in self.G.nodes:
            if "Latitude" in self.G.nodes[node_id].keys() and "hyperedge" in self.G.nodes[node_id].keys():
                latitudes.append(self.G.nodes[node_id]["Latitude"])
        return latitudes

    def get_latitude_limits(self):
        if self.latitude_limits == (-360, 360):
            latitudes = self.get_latitude_list()
            self.latitude_limits = (min(latitudes), max(latitudes))
            for edge in self.G.edges:
                _, edge_lats = self.edge_to_geo_proj(edge)
                if min(edge_lats) < self.latitude_limits[0]:
                    self.latitude_limits = (min(edge_lats), self.latitude_limits[1])
                if max(edge_lats) > self.latitude_limits[1]:
                    self.latitude_limits = (self.latitude_limits[0], max(edge_lats))
        return self.latitude_limits

    def get_longitude_list(self):
        longitudes = []
        for node_id in self.G.nodes:
            if "Longitude" in self.G.nodes[node_id].keys():
                longitudes.append(self.G.nodes[node_id]["Longitude"])
        return longitudes
    
    def get_hyperedge_longitude_list(self):
        longitudes = []
        for node_id in self.G.nodes:
            if "Longitude" in self.G.nodes[node_id].keys() and "hyperedge" in self.G.nodes[node_id].keys():
                longitudes.append(self.G.nodes[node_id]["Longitude"])
        return longitudes

    def get_longitude_limits(self):
        if self.longitude_limits == (-360, 360):
            logitudes = self.get_longitude_list()
            self.longitude_limits = (min(logitudes), max(logitudes))
            for edge in self.G.edges:
                edge_lons, _ = self.edge_to_geo_proj(edge)
                if min(edge_lons) < self.longitude_limits[0]:
                    self.longitude_limits = (min(edge_lons), self.longitude_limits[1])
                if max(edge_lons) > self.longitude_limits[1]:
                    self.longitude_limits = (self.longitude_limits[0], max(edge_lons))
        return self.longitude_limits
    
    def __calculate_hyperedge_locations(self):
        for node_id in self.G.nodes:
            node : dict = self.G.nodes[node_id]
            if "hyperedge" in node.keys():
                latitude = 0
                longitude = 0
                connection_count = 0
                for edge in self.G.edges:
                    connecting_node_id = -1
                    if node_id == edge[0]:
                        connecting_node_id = edge[1]
                    if node_id == edge[1]:
                        connecting_node_id = edge[0]
                    if connecting_node_id != -1:
                        if "Latitude" in self.G.nodes[connecting_node_id].keys():
                            connection_count += 1
                            latitude += self.G.nodes[connecting_node_id]["Latitude"]
                            longitude += self.G.nodes[connecting_node_id]["Longitude"]
                if connection_count == 0:
                    print(node_id)
                self.G.nodes[node_id]["Latitude"] = latitude / connection_count
                self.G.nodes[node_id]["Longitude"] = longitude / connection_count

    def plot_on_map(self, resolution = 'l'):
        self.__create_basemap_m(resolution)

        for node_id in self.get_nodes():
            self.__draw_node(self.get_node(node_id), refresh=False)


        for edge in self.G.edges:
            self.__draw_edge_on_geo(edge, refresh=False)

            


    def __draw_edge_on_geo(self, edge, color="green", refresh=True, inner_edge=False, in_italy = False):
            if edge in self.G.edges:
                if "points" in self.G.get_edge_data(edge[0], edge[1])[edge[2]].keys():
                    if not inner_edge:
                        # print(self.G.get_edge_data(edge[0], edge[1])[0]["points"]["point"])
                        for point_i, point in enumerate(self.G.get_edge_data(edge[0], edge[1])[edge[2]]["points"]["point"]):
                            if point_i > 0:
                                next_point = self.G.get_edge_data(edge[0], edge[1])[edge[2]]["points"]["point"][point_i - 1]
                                self.G.add_node("tmp1", Longitude=point["Longitude"], Latitude=point["Latitude"])
                                self.G.add_node("tmp2", Longitude=next_point["Longitude"], Latitude=next_point["Latitude"])
                                # self.__draw_node(self.get_node("tmp"))
                                # print("node drawn")
                                self.__draw_edge_on_geo(("tmp1", "tmp2"), inner_edge=True, refresh=False, in_italy=True)
                                self.G.remove_node("tmp1") 
                                self.G.remove_node("tmp2")    
            if "italy" not in self.file_name or in_italy:
                lons, lats = self.edge_to_geo_proj(edge, granularity=400)
                inflexion_i = 0
                xs, ys = [], []
                for i in range(len(lons)):
                    if i > 1 and ((lons[i] > 170 and lons[i-1] < -170) or (lons[i] < -170 and lons[i-1] > 170)):
                        inflexion_i = i
                        break
                    x, y = self.m(lons[i], lats[i])
                    xs.append(x)
                    ys.append(y)
                    
                plt.plot(xs, ys, color=color, linewidth=self.line_size * 2.5)
                if inflexion_i > 0:
                    xs, ys = [], []
                    for i in range(inflexion_i, len(lons)):
                        x, y = self.m(lons[i], lats[i])
                        xs.append(x)
                        ys.append(y)
                    plt.plot(xs, ys, color=color, linewidth=self.line_size * 2.5)
                if refresh:
                    plt.figure(1).canvas.draw()

    def __create_basemap_m(self, resolution):
        if "Continent" in self.G.graph["GeoExtent"]:
            if "Global" not in self.G.graph["Network"]:
                self.node_size = 1 # 20 for non large networks, .05 for global, 1 for large
                self.line_size = 0.1 # 1 for non large networks, .01 for global, .1 for large
            else:
                self.node_size = .01 # 20 for non large networks, .05 for global, 1 for large
                self.line_size = 0.05 # 1 for non large networks, .01 for global, .1 for large
        else:
            self.node_size = 10 # 20 for non large networks, .05 for global, 1 for large
            self.line_size = 0.3 # 1 for non large networks, .01 for global, .1 for large
        lat_limits = self.get_latitude_limits()
        lon_limits = self.get_longitude_limits()
        lat_diff = abs(lat_limits[0] - lat_limits[1])
        lon_diff = abs(lon_limits[0] - lon_limits[1])

        self.m = Basemap(
                projection='merc',  
                llcrnrlat=max(-88, lat_limits[0] - lat_diff * 0.2), urcrnrlat=min(88, lat_limits[1] + lat_diff * 0.2), 
                llcrnrlon=max(-180, lon_limits[0] - lon_diff * 0.2), urcrnrlon=min(180, lon_limits[1] + lon_diff * 0.2), 
                resolution=resolution
        )
        # self.m.etopo()
        # self.m.bluemarble()
        self.m.fillcontinents(color="lightgrey", lake_color="lightblue")
        self.m.drawmapboundary(color="white", fill_color="lightblue")
        self.m.drawcoastlines(linewidth=self.line_size*2)
        self.m.drawcountries(linewidth=self.line_size*2, color="white")  
        self.m.drawstates(linewidth=self.line_size, color="white")

    def edge_to_geo_proj(self, edge, granularity = 100):
        lon1, lat1 = self.G.nodes[edge[0]]["Longitude"], self.G.nodes[edge[0]]["Latitude"]
        lon2, lat2 = self.G.nodes[edge[1]]["Longitude"], self.G.nodes[edge[1]]["Latitude"]
        all_points = [(lon1, lat1)] + self.geod.npts(lon1, lat1, lon2, lat2, npts=granularity) + [(lon2, lat2)]
        lons, lats = zip(*all_points)
        lons = list(lons)
        lats = list(lats)
        lats.insert(0, lat1)
        lons.insert(0, lon1)
        lats.append(lat2)
        lons.append(lon2)

        return lons, lats

    def _prepare_atlantic_traversal_if_necessary(self):
        self.should_route = False
        for edge in self.G.edges:
            if self.G.nodes[edge[0]]["Longitude"] > -10 and self.G.nodes[edge[1]]["Longitude"] < -30:
                self.should_route = True
            elif self.G.nodes[edge[0]]["Longitude"] < -30 and self.G.nodes[edge[1]]["Longitude"] > -10:
                self.should_route = True
        if self.should_route:
            atlantic_nodes = [
                ("Offshore-NYC", {"Longitude": -72.5, "Latitude": 40.5, "internal": 1, "label": "Offshore-NYC", "router" : 1}),
                ("Hallifax", {"Longitude": -63.57, "Latitude": 44.65, "internal": 1, "label": "Hallifax", "router" : 1}),
                ("Mid-Atlantic", {"Longitude": -40, "Latitude": 50, "internal": 1,"label": "Mid-Atlantic", "router" : 1}),
                ("Offshore-Ireland", {"Longitude": -12, "Latitude": 52, "internal": 1, "label": "Offshore-Ireland", "router" : 1}),
                ("UK-Landing", {"Longitude": -3, "Latitude": 51.3, "internal": 1, "label": "UK-Landing", "router" : 1})
            
            ]
            self.G.add_nodes_from(atlantic_nodes)
            new_edges = []
            for edge in self.G.edges:
                
                if self.G.nodes[edge[0]]["Longitude"] > -10 and self.G.nodes[edge[1]]["Longitude"] < -30:
                    new_edges.append((edge[0], "UK-Landing"))
                    new_edges.append((edge[1], "Offshore-NYC"))
                elif self.G.nodes[edge[0]]["Longitude"] < -30 and self.G.nodes[edge[1]]["Longitude"] > -10:
                    new_edges.append((edge[1], "UK-Landing"))
                    new_edges.append((edge[0], "Offshore-NYC"))
                else:
                    new_edges.append(edge)
            
            self.G.remove_edges_from(list(self.G.edges))
            self.G.add_edges_from(new_edges)
            self.G.add_edge("Hallifax", "Offshore-NYC")
            self.G.add_edge("Hallifax", "Mid-Atlantic")
            self.G.add_edge("Offshore-Ireland", "Mid-Atlantic")
            self.G.add_edge("Offshore-Ireland", "UK-Landing")

    def shortest_distance_between_nodes(self):
        min_dist = 1000000000000000
        for node_label_1 in self.G.nodes:
            lon1 = self.G.nodes[node_label_1]["Longitude"]
            lat1 = self.G.nodes[node_label_1]["Latitude"]
            for node_label_2 in self.G.nodes:
                if node_label_1 != node_label_2:
                    lon2 = self.G.nodes[node_label_2]["Longitude"]
                    lat2 = self.G.nodes[node_label_2]["Latitude"]
                    _, _, distance = self.geod.inv(lon1, lat1, lon2, lat2)
                    distance /= 1000 # convert to km
                    if distance < min_dist:
                        min_dist = distance
        return min_dist

    def diameter_of_map_in_kms(self):
        lat_min, lat_max = self.get_latitude_limits()
        lon_min, lon_max = self.get_longitude_limits()
        _, _, distance = self.geod.inv(lon_min, lat_min, lon_max, lat_max)
        distance /= 1000 # convert to km
        return distance

    def save_map_plot(self, format='svg'):
        plt.figure(1, dpi=300)
        self.plot_on_map()
        name = self.G.graph["Network"].replace(" ", "_")
        plt.savefig("../networks/plots/"+ name + "." + format, format=format, bbox_inches="tight")
        plt.close(1)

    def show_map_plot(self):
        fig = plt.figure(1, dpi=300)
        self.plot_on_map()
        plt.show()
        plt.close(1)

    def draw_on_map_plot(self):
        fig = plt.figure(1, dpi=300)
        self.plot_on_map()
        _ = fig.canvas.mpl_connect('button_press_event', self.__onclick)
        _ = fig.canvas.mpl_connect('key_press_event', self.__on_key)
        plt.show()
        plt.close(1)

    def __on_key(self, event):
        if event.key == "G":
            self.serialize_to_gml()
        if event.key == "E":
            self.__add_drawing_edge()
        if event.key == "D":
            self.__delete_drawing_edge()

    def __delete_drawing_edge(self):
        edge = self.drawing_selected_edge
        if edge != None:
            self.__draw_edge_on_geo(edge, "red")
            self.G.remove_edge(edge[0], edge[1])
            self.__remove_line_of_edge(self.drawing_selected_edge)

            self.drawing_selected_edge = None

    def __remove_line_of_edge(self, edge):
        axs = plt.gca()
        for line in axs.lines:
            lon, lat = self.m(line.get_xdata()[20], line.get_ydata()[20], inverse=True)
            if self.__distance_between_point_and_edge(lon, lat, edge) < 1:
                line.remove()
        plt.figure(1).canvas.draw()

    def __add_drawing_edge(self):
        if len(self.drawing_selected_nodes) == 2:
            for node in self.drawing_selected_nodes:
                self.__draw_node(node)
            print(self.G.has_edge(self.get_node_id(self.drawing_selected_nodes[0]), self.get_node_id(self.drawing_selected_nodes[1])))
            if not self.G.has_edge(self.get_node_id(self.drawing_selected_nodes[0]), self.get_node_id(self.drawing_selected_nodes[1])):
                self.G.add_edge(self.get_node_id(self.drawing_selected_nodes[0]), self.get_node_id(self.drawing_selected_nodes[1]))
                self.__draw_edge_on_geo((self.get_node_id(self.drawing_selected_nodes[0]), self.get_node_id(self.drawing_selected_nodes[1])))
                plt.figure(1).canvas.draw()
            self.drawing_selected_nodes = []

    def get_node_id(self, node):
        for node_id in self.get_nodes():
            if self.get_node(node_id) == node:
                return node_id

    def __onclick(self, event):
        self.modified_by_drawing = True
        if event.inaxes:  # Only react if click happens inside the axes
            lon, lat = self.m(event.xdata, event.ydata, inverse=True)
            print(f"Clicked at data coords: ({lon}, {lat})")
            closest_node_dist = 100
            closest_edge_dist = 100
            closest_node = None
            closest_edge = None
            for node_id in self.get_nodes():
                node = self.get_node(node_id)
                _, _, distance = self.geod.inv(lon, lat, node["Longitude"], node["Latitude"])
                distance /= 1000 # convert to km
                if distance < closest_node_dist:
                    closest_node_dist = distance
                    closest_node = node
            for edge in self.G.edges:
                distance = self.__distance_between_point_and_edge(lon, lat, edge)
                if distance < closest_edge_dist:
                    closest_edge_dist = distance
                    closest_edge = edge

            if closest_edge_dist + 2 < closest_node_dist:
                closest_node = None

            if closest_node != None:
                self.drawing_selected_nodes.append(closest_node)
                self.__draw_node(closest_node, "green")
                if self.drawing_selected_edge != None:
                    self.__draw_edge_on_geo(self.drawing_selected_edge)
                    self.drawing_selected_edge = None
                if len(self.drawing_selected_nodes) > 2:
                    oldest_node = self.drawing_selected_nodes[0].copy()
                    self.drawing_selected_nodes = self.drawing_selected_nodes[1:]
                    print(oldest_node)
                    self.__draw_node(oldest_node)
            elif closest_edge != None:
                for node in self.drawing_selected_nodes:
                    self.__draw_node(node)
                if self.drawing_selected_edge != None:
                    self.__draw_edge_on_geo(self.drawing_selected_edge)
                    self.drawing_selected_edge = None
                self.drawing_selected_nodes = []
                self.__draw_edge_on_geo(closest_edge, color = "yellow")
                self.drawing_selected_edge = closest_edge

    def __distance_between_point_and_edge(self, lon, lat, edge):
        lons, lats = self.edge_to_geo_proj(edge)
        min_dist = 1000
        for i in range(len(lons)):
            _, _, distance = self.geod.inv(lon, lat, lons[i], lats[i])
            distance /= 1000
            if distance < min_dist:
                min_dist = distance
        return min_dist
                
    def __draw_node(self, node, color="red", refresh=True):
        lat = node["Latitude"]
        lon = node["Longitude"]
        if "hyperedge" in node.keys() and color == "red":
            color = "yellow"
        if "router" in node.keys() and color == "red":
            color = "purple"
        x, y = self.m(lon, lat)
        self.m.scatter(x, y, marker='.', color=color, zorder=6, s=self.node_size*2)
        if refresh:
            plt.figure(1).canvas.draw()

