import os
import shutil
import geopandas as gpd
import numpy as np
from scipy.spatial import KDTree
from shapely.geometry import *
import re
import json

# Define the directory where your GeoJSON files are located
walk_directory = '/Users/nathan/desktop/the_walk/data-python/walk'
gdf2_directory = '/Users/nathan/desktop/the_walk/data-python/filtered'
output_directory = '/Users/nathan/desktop/the_walk/data-python/output_knearest_geojson'  # New directory for output JSON file

# Delete the contents of the output directory if it exists
if os.path.exists(output_directory):
    shutil.rmtree(output_directory)

# Now create the output directory again
os.makedirs(output_directory)

# Continue with the rest of your script...

# Define a function to create a GeoJSON polygon feature from a list of coordinates
def create_polygon_feature(coords):
    return {
        "type": "Feature",
        "properties": {},
        "geometry": {
            "type": "Polygon",
            "coordinates": [coords]
        }
    }

# Get a list of all GeoJSON files in the directory
geojson_files = [f for f in os.listdir(walk_directory) if f.endswith('.geojson')]

# Load the first GeoJSON file into a GeoDataFrame
if geojson_files:
    gdf1 = gpd.read_file(os.path.join(walk_directory, geojson_files[0]))
else:
    raise FileNotFoundError("No GeoJSON files found in the specified directory.")

# Extract all points from LineString geometries in gdf1
walk_geom = np.concatenate([np.array(geom.coords) for geom in gdf1.geometry if isinstance(geom, LineString)])

# Get a list of all GeoJSON files in the directory
gdf2_geojson_files = [f for f in os.listdir(gdf2_directory) if f.endswith('.geojson')]

# Define a function to extract the desired names from the file paths
def get_name_from_file_path(file_path):
    file_name = os.path.splitext(os.path.basename(file_path))[0]  # Get the file name without extension
    if file_name.startswith("filtered_"):
        return file_name[len("filtered_"):]
    return file_name

# Create the gdf2_files list with dynamically constructed file paths and names
gdf2_files = [(os.path.join(gdf2_directory, file), get_name_from_file_path(file)) for file in gdf2_geojson_files]

# Load GeoJSON files into GeoDataFrames and extract coordinates with corresponding file names
gdf2_coords_with_names = []
for file_path, file_name in gdf2_files:
    gdf2 = gpd.read_file(file_path)
    gdf2_coords_with_names.extend([(file_name, coord) for geom in gdf2.geometry if isinstance(geom, Point) for coord in geom.coords])

# Extract coordinates and file names separately from the combined list
file_names, gdf2_coords = zip(*gdf2_coords_with_names)
gdf2_coords = np.array(gdf2_coords)
file_names = np.array(file_names)

# Set printing options to display full precision without truncation
np.set_printoptions(formatter={'float': '{:0.10f}'.format})

# Perform k-d tree nearest neighbor search for each point in walk_geom
tree = KDTree(gdf2_coords)

nearest_neighbors = []
for point in walk_geom:
    # Find fifteen nearest neighbors
    distances, indices = tree.query(point.reshape(1, -1), k=15)  # Set k to 15 for fifteen nearest neighbors
    # Get the corresponding file names and coordinates for the nearest neighbors using indices
    nearest_files = file_names[indices.flatten()]
    nearest_coords = gdf2_coords[indices.flatten()]
    nearest_neighbors.append((point, zip(nearest_files, nearest_coords)))

# Format the nearest neighbors data
features = []
for point, nearest_neighbors_data in nearest_neighbors:
    polygon_coords = [point.tolist()] + [coord.tolist() for _, coord in nearest_neighbors_data]
    feature = create_polygon_feature(polygon_coords)
    features.append(feature)

output_geojson = {
    "type": "FeatureCollection",
    "features": features
}

geojson_file_path = os.path.join(output_directory, 'output.geojson')
with open(geojson_file_path, 'w') as geojson_file:
    json.dump(output_geojson, geojson_file, indent=2)

print(f"Output GeoJSON file saved at: {geojson_file_path}")
