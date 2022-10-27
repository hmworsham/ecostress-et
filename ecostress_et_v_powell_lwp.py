# Import packages
import datetime as dt
import geopandas as gpd
import json
import matplotlib as mp
import numpy as np
import os
import pandas as pd
import pyproj
import rasterio as rio
import seaborn as sns

from ecostress_et import et

from matplotlib import pyplot as plt
from matplotlib import colors as colors
from matplotlib import cm as cmx

# Specify workspace
workdir = '/Volumes/GoogleDrive/My Drive/Research/ECOSTRESS/Data/'
indir = os.path.join(workdir, 'ET_PT-JPL_Point_Samples', 'EastRiver_2018-2021')
lwpdir = os.path.join(workdir, 'Powell_Eastriver_LWP')
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

# Prepare date and time values
dates = pd.to_datetime(eco['Date'])
eco['Date'] = dates.dt.date

# Create dataframe of leaf water potential data
lwp = pd.read_csv(os.path.join(
    lwpdir,
    [i for i in os.listdir(lwpdir) if 'csv' in i][0]),
    parse_dates={'Date': ['Year', 'Month', 'Day']})

# Summarize by date, tree, location, and species
lwpmean = lwp.groupby(['Date'], as_index=False).mean()

max(lwpmean.Date)

# Filter ec and eco by dates
eco1819 = eco[(eco['Date'] >= dt.date(2019, 6, 15)) & (eco['Date'] < dt.date(2019, 12, 15))]

eco1819 = eco1819[eco1819['ID']=='ER-APL1']

# Plot
# x = pd.to_datetime(eco1819['Date'])
# x = x.dt.strftime('%y-%m-%d')
# x = mp.dates.datestr2num(x)
# x = [dt.datetime.strptime(d, '%y-%m-%d') for d in eco1819['Date'].dt]

fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = plt.twinx()
sns.scatterplot(x=eco1819['Date'], y=eco1819['ECO3ETPTJPL_001_EVAPOTRANSPIRATION_PT_JPL_ETdaily'], label='ecostress instantaneous', marker='+')
sns.lineplot(x=lwpmean.Date, y=lwpmean.LWP, ax=ax2, label='predawn LWP')
# ax.scatter(
#     eco1819['Date'], eco1819['ECO3ETPTJPL_001_EVAPOTRANSPIRATION_PT_JPL_ETdaily'], marker='+')
#ax.scatter(ec1819['DATE'], ec1819['ET_W/m2/day'], marker='+')
ax.figure.legend()
plt.xticks(rotation=45)
ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m'))
ax.xaxis.set_major_locator(mpl.dates.WeekdayLocator(interval=3))
plt.show()

eco['ECO3ETPTJPL_001_EVAPOTRANSPIRATION_PT_JPL_ETdaily']
ec['ET_mm/day']
