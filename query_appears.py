# Query AppEARS 
#%% 
import geopandas as gpd
import requests as r
import getpass
import pprint
import time
import os
import cgi
import json

#%%
# Set indir, change wd
indir = '/Volumes/GoogleDrive/My Drive/Research/ECOSTRESS/Data/'
os.chdir(indir)
api = 'https://lpdaacsvc.cr.usgs.gov/appeears/api/'
outdir = os.path.join(indir, 'L3_ET_PT-JPL')

#%%
user = getpass.getpass(prompt = 'Enter NASA Earthdata Username: ')

password = getpass.getpass(prompt = 'Enter NASA Earthdata Password: ')

#%%
token_response = r.post('{}login'.format(api), auth=(user, password)).json()
del user, password
token_response

#%%
product_response = r.get('{}product'.format(api)).json()
productnames = [p['ProductAndVersion'] for p in product_response]
products = {p['ProductAndVersion']: p for p in product_response}

#%%
for p in productnames:
    if 'ECO' in products[p]['ProductAndVersion']:
        pprint.pprint(products[p])

# %%

# Define ecostress et and geolocation products
eco = 'ECO3ETPTJPL.001'
geo = 'ECO1BGEO.001'
eco_et = products[eco]
eco_geo = products[geo]

# %%

# Request layers for ecostress ET PT-JPL
et_response = r.get('{}product/{}'.format(api, eco)).json()
et_layers = list(et_response.keys())

# %%
# View metadata for a given layer
et_response['EVAPOTRANSPIRATION_PT_JPL_ETinst']

# Request layers for ecostress geolocation
geo_response = r.get('{}product/{}'.format(api, geo)).json()
geo_layers = list(geo_response.keys())

# %%
# Create tupled list linking desired product with desired layers

layers = [
    (eco, 'EVAPOTRANSPIRATION_PT_JPL_ETinst'), 
    (eco, 'EVAPOTRANSPIRATION_PT_JPL_ETinstUncertainty'),
    (eco, 'EVAPOTRANSPIRATION_PT_JPL_ETdaily'),
    (eco, 'EVAPOTRANSPIRATION_PT_JPL_ETcanopy')
    ]

for l in geo_layers:
    layers.append((geo, l))

#%% 
# Take the tupled list of layers and create a list of dictionaries to store each layer/product combo for ingestion in data request

prodlayer_ls = []
for l in layers:
    prodlayer_ls.append({
        'layer': l[1],
        'product': l[0]
    })
prodlayer_ls

# %%
# Request data
token = token_response['token']
reqhead = {'Authorization': 'Bearer {}'.format(token)}

# %%
# Define AOI for data request
aoi = gpd.read_file('../../RMBL/RMBL-East River Watershed Forest Data/Data/Geospatial/EastRiver_HU10_Extent')

#aoi = aoi.to_json()
#aoi = json.loads(aoi)

# %%
# Compile json task object

# Prompt for user-defined name of task
def appears_request(
    task_type, layers, proj, 
    out_format, start_date, end_date, 
    recurring=False, year_range = [2020,2022]):
    
    task_name = input('Enter a Task Name: ')

    task = {
        'task_type': task_type,
        'task_name': task_name,
        'params': {
            'dates': [
                {
                    'startDate': start_date,
                    'endDate': end_date,
                }],
                'layers': layers,
                'output': {
                    'format': {
                        'type': out_format},
                        'projection': proj},
                'geo': aoi,
                }
        }
    
    return task
# %%
task = appears_request(
        task_type='area',
        layers=prodlayer_ls,
        proj='geographic',
        out_format='geotiff', 
        start_date='04-15-2021',
        end_date='09-15-2021',
        recurring=False, 
        year_range=[2020, 2022]
    )

task_response = r.post(
    '{}task'.format(api),
    json = task,
    headers = reqhead).json()

# %%
statusparams = {'limit': 2, 'pretty': True} # Limit API response to 2 most recent entries, return as pretty json

tasks_response = r.get('{}task'.format(api), params=statusparams, headers=reqhead).json() # Query task service, setting params and header 

tasks_response

# %%
# Check status
task_id = '36b9fe9c-daf7-433f-94af-dea6627efcf6'
status_response = r.get('{}status/{}'.format(api, task_id), headers = reqhead).json()
status_response

#%%
# Ping call service for status
starttime = time.time()
while r.get('{}task/{}'.format(api, task_id), headers=reqhead).json()['status'] != 'done':
    print(r.get('{}task/{}'.format(api, task_id), headers=reqhead).json()['status'])
    time.sleep(30.0 - ((time.time() - starttime) % 30.0))
print(r.get('{}task/{}'.format(api, task_id), headers=reqhead).json()['status'])

# %%
# examine files in bundle
bundle = r.get('{}bundle/{}'.format(api,task_id)).json()

# %%
files = {}
for f in bundle['files']: files[f['file_id']] = f['file_name']
files

# %%
def download_appeears(files):

    for f in files:
        
        # Get a stream to the bundle file
        dl = r.get('{}bundle/{}/{}'.format(api, task_id, f), stream=True)
        
        #Parse the name from Content-Disposition header 
        filename = os.path.basename(cgi.parse_header(dl.headers['Content-Disposition'])[1]['filename'])
        
        # Create output file path
        filepath = os.path.join(outdir, filename)
        with open(filepath, 'wb') as f:
            
            # Write file to dest dir
            for data in dl.iter_content(chunk_size=8192): f.write(data)
    
    print('Downloaded files can be found at: {}'.format(outdir))

download_appeears(files)
# %%
