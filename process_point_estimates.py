# Import packages
import geopandas as gpd
import json
import matplotlib as mpl
import numpy as np
import os
import pandas as pd
import pyproj
import rasterio as rio

from ecostress_et import et

from matplotlib import pyplot as plt
from matplotlib import colors as colors
from matplotlib import cm as cmx

# Specify workspace
workdir = '/Volumes/GoogleDrive/My Drive/Research/ECOSTRESS/Data/'
indir = os.path.join(workdir, 'ET_PT-JPL_Point_Samples', 'EastRiver_2018-2021')
outdir = os.path.join(workdir, 'Output')

# List files in results of ecostress request
ecofiles = os.listdir(indir)
ecocsv = [i for i in ecofiles if 'csv' in i]
ecotxt = [i for i in ecofiles if 'txt' in i]
ecojson = [i for i in ecofiles if 'json' in i]

# Create list of csv files
ecodfs = []
for i in range(len(ecocsv)):
    f = pd.read_csv(os.path.join(indir, ecocsv[i]))
    ecodfs.append(f)

# Create dataframe of all ecostress results
i = 0
eco = ecodfs[i]
while i < 3:
    eco = eco.merge(
        ecodfs[i+1],
        how='left',
        on=['ID', 'Date'],
        suffixes=('', '_DROP')).filter(regex='^(?!.*_DROP)')
    i += 1

cols = '\n'.join(str(i) for i in eco.columns)

with open('/Users/hmworsham/Desktop/ecostress_columns.txt', 'w') as f:
    for lines in cols:
        f.write(lines)

# Create list of granules
with open(os.path.join(indir, ecotxt[0])) as g:
    granules = g.read()
granules = granules.split('\n')

# Create dataframe of points of observation
with open(os.path.join(indir, ecojson[0])) as j:
    req = json.load(j)

points = pd.DataFrame(req['params']['coordinates'])

##################################
# Filter observations by QA flags
##################################
vz = 'ECO1BGEO_001_Geolocation_view_zenith'
va = 'ECO1BGEO_001_Geolocation_view_azimuth'
aq = 'ECO3ANCQA_001_L3_L4_QA_MOD13Q1_QC_Aerosol_Quantity'
gh = 'ECO1BGEO_001_Geolocation_height'
uc = 'ECO3ETPTJPL_001_EVAPOTRANSPIRATION_PT_JPL_ETinstUncertainty'
ee = eco.dropna(0, subset=[vz, uc])

xx = ee[vz]
yy = ee[uc]
plt.scatter(xx,yy)
from sklearn.linear_model import LinearRegression
XX = xx.values.reshape(-1,1)
YY = yy.values.reshape(-1,1)
reg = LinearRegression().fit(XX, YY)
slope = reg.coef_
r2 = reg.score(XX,YY)
r2

eco = eco[eco[vz]<=20]

#######################
# Plots
#######################
# Plot granules by date and overpass time in Mountain Time
# Prepare date and time values
dates, times = et.set_datetime(eco, tz='US/Mountain', unique=True)
dates.shape
x, y = et.datetime2num(dates, times)
eco['Date'] = x
eco['Time'] = y

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.axvspan(dt.datetime(2018, 4, 1), dt.datetime(2019, 1, 8), color='orange', alpha=0.9)
ax.axvspan(dt.datetime(2020, 6, 2), dt.datetime(2021, 7, 21), color='orange', alpha=0.9)
ax.scatter(eco.Date, eco.Time, marker='+')
plt.xticks(rotation=45)
ax.yaxis.set_major_locator(mpl.dates.HourLocator(interval=2))
ax.yaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
ax.set_xlabel('Overpass Date')
ax.set_ylabel('Overpass Time (Zone: US-Mountain)')

# Set the color map to match the number of species
ids = pd.unique(eco.ID)
z = range(1, len(ids))
hot = plt.get_cmap('viridis')
cnorm = colors.Normalize(vmin=0, vmax=len(ids))
scalarmap = cmx.ScalarMappable(norm=cnorm, cmap=hot)
scalarmap.to_rgba(1)
