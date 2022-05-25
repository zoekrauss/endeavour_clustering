from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from eqcorrscan.utils.catalog_utils import filter_picks
import glob2 as glob
import numpy as np
import obspy
import matplotlib.pyplot as plt
from eqcorrscan.core.template_gen import template_gen


## READ IN CATALOG ##
# Read in quakeML file
from obspy.core.event import read_events
cat = read_events('endquakes_july2019.xml')
cat[0].origins[0].arrivals
pick_id = cat[0].picks[0].resource_id
arrivals = cat[0].origins[0].arrivals
arr = [a for a in arrivals if a.pick_id==pick_id]
# Add in a phase hint to the picks
for event in cat.events:
    for pick in event.picks:
        pick_id = pick.resource_id
        arr = [a for a in event.origins[0].arrivals if a.pick_id==pick_id]
        pick.phase_hint=arr[0].phase
print('Read the catalog!')
        
        
## FILTER CATALOG ##
# Filter the catalog to only vent field earthquakes, for which the window length 
# should be 2-2.5 s judging by visual examination of the waveforms
latlon = [47.87,48.05,-129.25,-128.97]
cat_filt = obspy.core.event.Catalog()
for ev in cat:
    lat = ev.origins[0].latitude
    lon = ev.origins[0].longitude
    if (lat>latlon[0]) & (lat < latlon[1]) & (lon>latlon[2]) & (lon<latlon[3]):
        cat_filt.append(ev)
        

## READ IN TEMPLATES FROM IRIS ##
templates=template_gen(method="from_client",catalog=cat_filt,client_id='IRIS',lowcut=5,highcut=20,samp_rate=200.0, filt_order=4, length=2.0,prepick=0.1, swin='all',process_len=200,parallel=True)
template_list = []
for i,t in enumerate(templates):
    template_list.append([t,cat_filt[i].resource_id.id])
print('Read in templates from IRIS!')

## PERFORM THE CROSS CORRELATION OF THE TEMPLATES ##
from obspy import read
import glob
import os
from eqcorrscan.utils.clustering import cluster
from eqcorrscan import tests
# Compute cross-correlation coefficients between all templates
import eqcorrscan
dist,shift,ids=eqcorrscan.utils.clustering.distance_matrix(templates,shift_len=0.0,allow_individual_trace_shifts=True, cores=12)
print('Done with cross correlation, time to save!')
# Save the results
np.save('july2019_dist',dist)
np.save('july2019_shift',shift)
np.save('july2019_ids',ids)