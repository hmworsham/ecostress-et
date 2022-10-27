# %%
# Import packages
import datetime as dt
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

# %%
# Specify workspace
workdir = '/Volumes/GoogleDrive/My Drive/Research/ECOSTRESS/Data/'

indir = os.path.join(
    workdir,
    'ET_PT-JPL_Point_Samples',
    'ryken-ec-ecostress_2018-2021')
ecdir = os.path.join(
    workdir,
    'Ryken_EastRiver_EC')
outdir = os.path.join(
    workdir,
    'Output')
# %%
# List files in results of ecostress request
ecofiles = os.listdir(indir)
ecocsv = [i for i in ecofiles if 'csv' in i]
ecotxt = [i for i in ecofiles if 'txt' in i]
ecojson = [i for i in ecofiles if 'json' in i]

# Create list of csv files
ecodfs = []
for i in range(len(ecocsv)):
    f = pd.read_csv(os.path.join(indir, ecocsv[i]),na_values=9999)
    ecodfs.append(f)

# Create dataframe of all ecostress results
i = 0
eco = ecodfs[i]
while i < 6:
    eco = eco.merge(
        ecodfs[i+1],
        how='left',
        on=['ID', 'Date'],
        suffixes=('', '_DROP')).filter(regex='^(?!.*_DROP)')
    i += 1

#%%
cols = '\n'.join(str(i) for i in eco.columns)

with open('/Users/hmworsham/Desktop/ecostress_columns.txt', 'w') as f:
    for lines in cols:
        f.write(lines)

# %%
# Prepare date and time values
dates = pd.to_datetime(eco['Date'])
eco['Date'] = dates.dt.date
eco['Date'] = pd.to_datetime(eco['Date'])

# %%
# Create dataframe of eddy covariance data
ec = pd.read_csv(os.path.join(ecdir, 'FluxTower_Pumphouse_ESS-DIVE.csv'))

ec['DATE'] = pd.to_datetime(ec.DAY_MT)
ec['ET_W/m2/day'] = ec['ET_mm/day']/0.0864/0.408

# %%
# Filter ec and eco by dates
ec1819 = ec[(ec['DATE'] >= dt.datetime(2018, 4, 1)) &
            (ec['DATE'] < dt.datetime(2020, 1, 1))]
eco1819 = eco[(eco['Date'] >= dt.datetime(2018, 4, 1)) &
              (eco['Date'] < dt.datetime(2020, 1, 1))]

# ec1819 = ec[(ec['DATE'] >= min(eco['DATE'])) &
#             (ec['DATE'] < max(ec['DATE']))]
# eco1819 = eco[(eco['Date'] >= min(eco['DATE']) &
#               (eco['Date'] < min(eco['DATE'])]

# %%
# Set variable shortcuts
et_inst = 'ECO3ETPTJPL_001_EVAPOTRANSPIRATION_PT_JPL_ETinst'
et_daily = 'ECO3ETPTJPL_001_EVAPOTRANSPIRATION_PT_JPL_ETdaily'
eco1819[et_inst]
et_alexi_daily = 'ECO3ETALEXI_001_EVAPOTRANSPIRATION_ALEXI_ETdaily'


#%%
eco[eco1819[et_alexi_daily]9999]

# %%
# Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(
    ec1819['DATE'], 
    ec1819['LE_W/m2'], 
    c='lightgreen', 
    s=2,
    label = 'Eddy Covariance Daily ET')
ax.scatter(
    eco1819['Date'], 
    eco1819[et_daily], 
    c='dodgerblue', 
    marker='+', 
    label='ECOSTRESS PT-JPL Daily ET')
ax.scatter(
    eco1819['Date'], 
    eco1819[et_alexi_daily], 
    c='navy', 
    marker='+',
    label='ECOSTRESS disALEXI Daily ET')
plt.xticks(rotation=45)
ax.legend()
ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m'))
ax.xaxis.set_major_locator(mpl.dates.WeekdayLocator(interval=6))
ax.set_xlabel('Observation Date')
ax.set_ylabel('Evapotranspiration (W/m2)')


# %%
# Filter ec and eco by dates
# ecset = ec[ec['DATE'].isin(list(eco['Date']))]
# ecoset = eco[eco['Date'].isin(list(ecset['DATE']))]

ecoec = ec.merge(eco, left_on='DATE', right_on='Date')
ecoec = ecoec.dropna(0, subset=['LE_W/m2', et_inst])
# %%
x=ecoec['LE_W/m2']
y=ecoec[et_daily]

# %%
X = x.values.reshape(-1,1)
Y = y.values.reshape(-1,1)
reg = LinearRegression().fit(X, Y)
slope = reg.coef_
r2 = reg.score(X,Y)

#%%
fig = plt.figure(figsize=((8,8)))
#ax = fig.add_subplot(111)
ax = plt.gca()
plt.scatter(x,y, c='dodgerblue', marker='+', s=40)
plt.plot(np.arange(300), np.arange(300), c='k', linewidth=0.5)
plt.xlim([0,300])
plt.ylim([0,300])
plt.text(.75,.88,
        'Slope={:.2g}\nR²{:.2g}'.format(float(slope),float(r2)),
        transform=ax.transAxes,
        fontsize=12,
        ha='left')
plt.show()
# %%
from sklearn.linear_model import LinearRegression

# %%
