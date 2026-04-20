from pyproj import Geod
import geopandas as gpd
from shapely.geometry import Point, LineString

# naturalearthdata.com -> "Physical" -> "Land" 1:10m
land = gpd.read_file("../inputs/ne_10m/ne_10m_land.shp")  # path to shapefile


geod = Geod(ellps="WGS84")

def sample_geodesic(lon1, lat1, lon2, lat2, spacing_m=1000):
    # compute total geodesic distance and initial azimuth
    az12, az21, total_m = geod.inv(lon1, lat1, lon2, lat2)
    # how many segments
    if total_m == 0:
        return [(lon1, lat1)]
    n = max(1, int(total_m // spacing_m))
    pts = geod.npts(lon1, lat1, lon2, lat2, n-1) if n > 1 else []
    pts = [(lon1, lat1)] + pts + [(lon2, lat2)]
    return pts

# Test whether the geodesic crosses water
def crosses_water(lon1, lat1, lon2, lat2, spacing_m=1000, land_gdf=land):
    pts = sample_geodesic(lon1, lat1, lon2, lat2, spacing_m)
    land_sindex = land_gdf.sindex
    for lon, lat in pts:
        p = Point(lon, lat)
        possible_matches_index = list(land_sindex.intersection(p.bounds))
        if not possible_matches_index:
            return True
        if not any(land_gdf.iloc[i].geometry.contains(p) for i in possible_matches_index):
            return True
    return False



