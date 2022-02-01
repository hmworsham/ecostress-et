# Import libraries
import os
import geopandas as gpd
import numpy as np
import pandas as pd
import datetime
import rasterio as rio

from osgeo import gdal, osr, gdalconst
from rasterio import mask
from matplotlib import pyplot as plt

# Set up workspace
datadir = '/Volumes/GoogleDrive/My Drive/Research/'
ecodir = os.path.join(datadir, 'ECOSTRESS/Data/Output')
topodir = os.path.join(datadir, 'RMBL/RMBL-East River Watershed Forest Data/Data/Geospatial/Worsham_2021_SiteSelection/2021_Analysis_Layers/USGS_1-9_arcsec_DEM')

# 
topopaths = []
for root,sub,files in os.walk(topodir):
    for f in files:
        if f.endswith('tif'):
            topopaths.append(os.path.join(root,f))

with rio.open(topopaths[2], 'r') as src:
    data = src.read(1)
    profile = src.profile
    tags = src.tags()

topos = []
profiles = []
tags = []

maskedarr = np.ma.array(topos[0], mask=(topos[0]<=0))
for tif in topopaths:
    with rio.open(topopaths[0], 'r') as src:
        data = src.read()
        profile = src.profile
        tag = src.tags()
        name = tif.split('/')[-1]
        topos.append(data)
        profiles.append(profile)
        tags.append(tag)

def cliprasters(raster):
    
    geompath = os.path.join(datadir, 'RMBL/RMBL-East River Watershed Forest Data/Data/Geospatial/EastRiver_HU10_Extent')

    aoi = gpd.read_file(geompath).geometry
    
    with rio.open(raster) as src:
        aoi = aoi.to_crs(src.profile['crs'])
        out_image, out_transform = rio.mask.mask(src, aoi, crop=True)
        out_image = np.ma.array(out_image, mask=out_image[0] <= 0)

        out_meta = src.meta

        
        
    out_meta.update({"driver": "GTiff",
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform})
    
    return out_image, out_meta



img1, meta1 = cliprasters(topopaths[0], geompath)

plt.imshow(img1[0], cmap='gray', vmin=0, vmax=255)


topo_clips = list(map(cliprasters, topopaths))

# List evapotranspiration files
etpaths = [os.path.join(ecodir, i) for i in os.listdir(ecodir) if 'xml' not in i and 'tif' in i]

# Define path to a specific file
et_clip = list(map(cliprasters, etpaths))


# Co-register
reference = et_clip[0]
referenceProj = reference[1]['crs']
referenceTrans = reference[1]['transform']
referenceDtype = reference[1]['dtype']
x = reference[1]['width']
y = reference[1]['height']

rio.warp.aligned_target(referenceTrans, x, y, 0.0006302015699066588)


destination = np.zeros(reference[0].shape, dtype=np.uint8)
dst_crs = "EPSG:4326"

_, dst_transform = rreproject(
    source,
    destination,
    rpcs=source.rpcs,
    src_crs=src_crs,
    dst_crs=dst_crs,
    resampling=Resampling.nearest,
    **kwargs
)

# reference = et_clip[0]
# referenceProj = reference[1]['crs']
# referenceTrans = reference[1]['transform']
# referenceDtype = reference[1]['dtype']
# x = reference[1]['width']
# y = reference[1]['height']

# srs = osr.SpatialReference()
# srs.ImportFromEPSG(4326)
# srs = srs.ExportToWkt()

# outputfile = os.path.join(ecodir, 'coreg1')
# driver= gdal.GetDriverByName('GTiff')
# output = driver.Create(outputfile, x, y, 1, gdal.GDT_Float32)
# output.SetGeoTransform(list(referenceTrans[:6]))
# output.SetProjection(srs)

# gdal.ReprojectImage(reference[0], output, referenceProj, referenceProj, gdalconst.GRA_Bilinear)


# OutTile = gdal.Warp(OutTileName, Raster,
#   format=RasterFormat, outputBounds=[minX, minY, maxX, maxY], 
#   xRes=PixelRes, yRes=PixelRes, dstSRS=Projection, 
#   resampleAlg=gdal.GRA_NearestNeighbour, options=['COMPRESS=DEFLATE']
# )
