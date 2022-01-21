# Import packages
import os
import pyproj
import geopandas as gpd
import numpy as np
import pandas as pd
import datetime
import rasterio as rio

from matplotlib import pyplot as plt

# Specify workspace
workdir = '/Volumes/GoogleDrive/My Drive/Research/ECOSTRESS/Data/'
indir = os.path.join(workdir, 'L3_ET_PT-JPL')
outdir = os.path.join(workdir, 'Output')

# Reorganize directory
filetypes = [
    'Geolocation_height',
    'Geolocation_land_fraction',
    'Geolocation_solar_azimuth',
    'Geolocation_solar_zenith',
    'Geolocation_view_azimuth',
    'Geolocation_view_zenith',
    'ETcanopy',
    'ETdaily',
    'ETinstUncertainty',
    'ETinst'
]

# Create subdirectories for each filetype and move files of each type to its new subdirectory.
for ft in filetypes:
    fp = os.path.join(outdir, ft)
    if not os.path.isdir(fp):
        print('Creating directory: {}'.format(fp))
        os.mkdir(fp)
        targets = [f for f in os.listdir(indir) if ft in f]
        inds = 0
        for i in range(len(targets)):
            oldpath = os.path.join(indir, targets[i])
            newpath = os.path.join(fp, targets[i])
            inds+=1
            os.rename(oldpath, newpath)
        print('{} files moved to {}'.format(inds,fp))
    else:
        print('Directory "{}" already exists'.format(fp))


# List evapotranspiration files
etfiles = os.listdir(os.path.join(indir, 'ETinst'))
etfiles = [i for i in etfiles if 'xml' not in i and 'tif' in i]

# Define path to a specific file
etpath = os.path.join(indir, 'ETinst', etfiles[18])

# Open file with rasterio
with rio.open(etpath, 'r') as src:
    et = src.read(1)
    profile = src.profile
    tags = src.tags()

# Pull units
try:
    units = tags['units']
except:
    units = 'none'

# Plot file
fig = plt.figure(figsize=(14,12))
im = plt.imshow(et, cmap='BuGn')
plt.colorbar(im, orientation='horizontal', fraction=0.05, pad=0.004, label=f'ET ({units})', shrink=0.6).outline.set_visible(True)
plt.tight_layout=True
plt.show()