#%%
# Import packages
import h5py
import os
from os.path import join
import pyproj
import geopandas as gpd
import numpy as np
import pandas as pd
import datetime
from dateutil import parser
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pyresample import geometry as geom
from pyresample import kd_tree as kdt
from osgeo import gdal, gdal_array, gdalconst, osr

# %%
# Set up working environment

# Current working directory will be set as the input directory
indir = '/Volumes/GoogleDrive/My Drive/Research/ECOSTRESS/Data/L3_ET_PT-JPL'
os.chdir(indir)
outdir = os.path.join('../Output')

# Create output directory
if not os.path.exists(outdir): 
    os.makedirs(outdir)

#%%
# List directory contents and create lists of ECOSTRESS HDF5 files (GEO, ET)
geo_list = [file for file in os.listdir(indir) if file.endswith('.h5') and 'GEO' in file]
print("geolocation:\n{}".format("\n".join(geo_list)))
eco_list = [file for file in os.listdir(indir) if file.endswith('.h5') and 'GEO' not in file]
print("products:\n{}".format("\n".join(eco_list)))

# %%
# Open HDF5 file
f = h5py.File(eco_list[0])             # Read in ECOSTRESS HDF5 file
eco_name = eco_list[0].split('.h5')[0]  # Keep original filename
print(eco_name)

# %%
# Create list of SDS inside HDF5 file
eco_objs = []
f.visit(eco_objs.append)
eco_sds = [str(obj) for obj in eco_objs if isinstance(f[obj], h5py.Dataset)] 

# Subset list to ETinst and ETinstUncertainty
sds = ['ETinst', 'ETinstUncertainty']
eco_sds = [dataset for dataset in eco_sds if dataset.endswith(tuple(sds))]
for dataset in eco_sds:
    print(dataset.split('/')[-1])
# %%
# import geolocation file
# Find the matching ECO1BGEO file from the file list
geo = [geo_file for geo_file in geo_list if eco_list[0][-37:-5] in geo_file]


# %%
# read in ECO1BGEO file, search for lat/lon SDS and import into Python as arrays
# Open Geo File
g = h5py.File(geo[0])
geo_objs = []
g.visit(geo_objs.append)

# Search for lat/lon SDS inside data file
lat_sd = [str(obj) for obj in geo_objs if isinstance(g[obj], h5py.Dataset) and '/latitude' in obj]
lon_sd = [str(obj) for obj in geo_objs if isinstance(g[obj], h5py.Dataset) and '/longitude' in obj]

# Open SDS as arrays
lat = g[lat_sd[0]][()].astype(float)
lon = g[lon_sd[0]][()].astype(float)

# Read the array dimensions
dims = lat.shape
print(dims)

# %%
# Set swath definition from lat/lon arrays
swath_def = geom.SwathDefinition(lons=lon, lats=lat)
swath_def.corners
# %%

# Define the lat/lon for the middle of the swath
mid = [int(lat.shape[1] / 2) - 1, int(lat.shape[0] / 2) - 1]
midLat, midLon = lat[mid[0]][mid[1]], lon[mid[0]][mid[1]]
midLat, midLon

# %%

'''
Perform a cartographic transformation by defining an Azimuthal Equidistant projection centered on the midpoint of the swath. Once the projection is defined, convert the lower left and upper right corners of the lat/lon arrays to a location (in meters) in the new projection. Lastly, measure the distance between the corners and divide by 70 (meters), the nominal pixel size that we are aiming for. Azimuthal Equidistant projection was chosen here based on the following characteristics of this projection:
# - Units in meters (necessary for defining 70 m pixels)  
# - Distances between all points are proportionally correct from center point
# - Azimuth (direction) are correct from the center point
'''

# Define AEQD projection centered at swath center
epsg_convert = pyproj.Proj("+proj=aeqd +lat_0={} +lon_0={}".format(midLat, midLon))

# Use info from AEQD projection bbox to calculate output cols/rows/pixel size
ll_lon, ll_lat = epsg_convert(np.min(lon), np.min(lat), inverse=False)
ur_lon, ur_lat = epsg_convert(np.max(lon), np.max(lat), inverse=False)
area_extent = (ll_lon, ll_lat, ur_lon, ur_lat)
cols = int(round((area_extent[2] - area_extent[0]) / 70))  # 70 m pixel size
rows = int(round((area_extent[3] - area_extent[1]) / 70))

# %%
# Use number of rows and columns generated above from the AEQD projection to set a representative number of rows and columns in the Geographic  area definition, which will then be translated to degrees below, then take the smaller of the two pixel dims to determine output size and ensure square pixels.

# Define Geographic projection
epsg, proj, p_name = '4326', 'longlat', 'Geographic'

# Define bounding box of swath
ll_lon, ll_lat, ur_lon, ur_lat = np.min(lon), np.min(lat), np.max(lon), np.max(lat)
area_extent = (ll_lon, ll_lat, ur_lon, ur_lat)

# Create area definition with estimated number of columns and rows
proj_dict = {'proj': proj, 'ellps': 'WGS84', 'datum': 'WGS84'}
area_def = geom.AreaDefinition(epsg, p_name, proj, proj_dict, cols, rows, area_extent)
area_def
# %%
# Square the pixels
# Square pixels and calculate output cols/rows
ps = np.min([area_def.pixel_size_x, area_def.pixel_size_y])
cols = int(round((area_extent[2] - area_extent[0]) / ps))
rows = int(round((area_extent[3] - area_extent[1]) / ps))

# Set up a new Geographic area definition with the refined cols/rows
area_def = geom.AreaDefinition(epsg, p_name, proj, proj_dict, cols, rows, area_extent)

#%%
# Get arrays with information about the nearest neighbor to each grid point 
'''
This is the most computationally heavy task in the swath2grid conversion and using `get_neighbour_info` speeds up the process if you plan to resample multiple SDS within an ECOSTRESS product (compute once instead of for every SDS). 
'''

index, outdex, index_arr, dist_arr = kdt.get_neighbour_info(swath_def, area_def, 210, neighbours=1)

# %%
# List the attributes for the `ETinst` layer, which can then be used to define the fill value and scale factor.

# Read in ETinst and print out SDS attributes
s = eco_sds[0]
ecoSD = f[s][()] 
for attr in f[s].attrs:
    if type(f[s].attrs[attr]) == np.ndarray:
        print(f'{attr} = {f[s].attrs[attr][0]}')
    else:
        print(f'{attr} = {f[s].attrs[attr].decode("utf-8")}')
# %%
# Extract scale factor, add offset and impose fill value from SDS metadata
# Read SDS attributes and define fill value, add offset, and scale factor if available
try:
    fv = int(f[s].attrs['_FillValue'])
except KeyError:
    fv = None
except ValueError:
    fv = f[s].attrs['_FillValue'][0]
try:
    sf = f[s].attrs['_Scale'][0]
except:
    sf = 1
try:
    add_off = f[s].attrs['_Offset'][0]
except:
    add_off = 0
try:
    units = f[s].attrs['units'].decode("utf-8")
except:
    units = 'none'

#%% 
# Perform K-D Tree nearest neighbor resampling (swath 2 grid conversion)
# Remember that the resampling has been split into two steps. In section 3b. arrays containing the nearest neighbor to each grid point were created. The second step is to use those arrays to retrieve a resampled result. 
et_geo = kdt.get_sample_from_neighbour_info('nn', area_def.shape, ecoSD, index, outdex, index_arr, fill_value=fv)

# Define geotransform
# Define the geotransform 
gt = [area_def.area_extent[0], ps, 0, area_def.area_extent[3], 0, -ps]
gt
# %%

# 3e. Basic image processing
# Apply Scale Factor and Add Offset
et_geo = et_geo * sf + add_off

# Set Fill Value   
et_geo[et_geo == fv * sf + add_off] = fv  

# %%

# Run steps 3c-3e for ET uncertainty
s = eco_sds[1]
ecoSD = f[s][()]
try:
    fv = int(f[s].attrs['_FillValue'])
except KeyError:
    fv = None
except ValueError:
    fv = f[s].attrs['_FillValue'][0]
try:
    sf = f[s].attrs['_Scale'][0]
except:
    sf = 1
try:
    add_off = f[s].attrs['_Offset'][0]
except:
    add_off = 0
unc_geo = kdt.get_sample_from_neighbour_info(
    'nn', area_def.shape, ecoSD, index, 
    outdex, index_arr, fill_value=fv)

unc_geo = unc_geo * sf + add_off
unc_geo[unc_geo == fv * sf + add_off] = fv

# %%

####################
# Exporting results
####################

# Set up dictionary of arrays to export
out_files = {'ET_inst': et_geo, 'ET_inst_uncertainty': unc_geo}

# %%
# Define CRS and export as geotiffs
# Loop through each item in dictionary created above
for file in out_files:
    
    # Set up output name
    out_name = join(outdir, '{}_{}.tif'.format(eco_name, file))
    print("output file:\n{}\n".format(out_name))
    
    # Get driver, specify dimensions, define and set output geotransform
    height, width = out_files[file].shape
    driv = gdal.GetDriverByName('GTiff')
    data_type = gdal_array.NumericTypeCodeToGDALTypeCode(out_files[file].dtype)
    d = driv.Create(out_name, width, height, 1, data_type)
    d.SetGeoTransform(gt)
        
    # Create and set output projection, write output array data
    # Define target SRS
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg))
    d.SetProjection(srs.ExportToWkt())
    srs.ExportToWkt()
    
    # Write array to band
    band = d.GetRasterBand(1)
    band.WriteArray(out_files[file])
    
    # Define fill value if it exists, if not, set to mask fill value
    if fv is not None and fv != 'NaN' and ~np.isnan(fv):
        print(fv)
        band.SetNoDataValue(fv)
    else:
        try:
            band.SetNoDataValue(out_files[file].fill_value)
        except AttributeError:
            pass
    band.FlushCache()
    d, band = None, None    
# %%
########
# Subset to local point (e.g. EC tower)
########
tower_lat, tower_lon = 38.92448618, -107.00236263  # AmeriFlux US-CZ3 tower location

# Calculate tower lat/lon distance from upper left corner, then divide by pixel size to find x,y pixel location
tcol = int(round((tower_lon - gt[0]) / gt[1]))
trow = int(round((tower_lat - gt[3]) / gt[5]))

# Print ET at the tower location
et_geo[trow, tcol]  

et_footprint = et_geo[(trow - 100):(trow + 200), (tcol - 100):(tcol + 200)] 
unc_footprint = unc_geo[(trow - 100):(trow + 200), (tcol - 100):(tcol + 200)] 



############
# Summary and plotting
############

# %%
et_median = np.nanmedian(et_footprint)
unc_median = np.nanmedian(unc_footprint)
print(f"Median ET: {et_median:0.3f} \nUncertainty: {unc_median:0.3f}")

#%%
# Calculate local overpass time
# Grab UTC time of observation from file name
eco_time = eco_name.split('_')[-3]
eco_time

# Convert UTC time to local overpass
observationTime = parser.parse(eco_time)
solarOverpass = observationTime + datetime.timedelta(hours=(np.radians(tower_lon) / np.pi * 12))
overpass = datetime.time(solarOverpass.hour, solarOverpass.minute)
date = observationTime.strftime('%Y-%m-%d')

# %%
# KDE

title = 'ECO3ETPTJPL Evapotranspiration'
sds_name = eco_sds[0].split("/")[-1]

pd.DataFrame({title: et_subset.flatten()}).plot.kde() 
plt.title(f'Probability Density of {sds_name} Surrounding US-CZ3\n{date} at {solarOverpass}');
plt.xlabel(f'{sds_name} ({units})')
plt.show()
# %%

# Set a Radius and calculate subset region from flux tower location (row, col)
radius = 300
et_subset = et_geo[(trow - radius):(trow + radius + 1), (tcol - radius):(tcol + radius + 1)]


# %%

# Create a colormap for the ET data
et_colors = ["#f6e8c3", "#d8b365", "#99974a", "#53792d", "#6bdfd2", "#1839c5"]
et_cmap = LinearSegmentedColormap.from_list("ET", et_colors)

# Plot
fig, ax = plt.subplots(figsize=(14,12))
fig.suptitle(f'{title} ({sds_name})\n{date} at {solarOverpass}', fontsize=26)
#plt.axis('off')
plt.imshow(et_subset, cmap=et_cmap)
plt.scatter(et_subset.shape[0]/2, et_subset.shape[1]/2, color="black", marker='x')

# Add a colormap legend
plt.colorbar(im, orientation='horizontal', fraction=0.05, pad=0.004, label=f'ET ({units})', shrink=0.6, ax=ax).outline.set_visible(True)

#%%
# Set up file name and export to png file
figure_filename = join(outdir, "{}_{}.png".format(eco_name, sds_name))
print("figure filename: {}".format(figure_filename))
fig.savefig(figure_filename, dpi=300)
plt.show()