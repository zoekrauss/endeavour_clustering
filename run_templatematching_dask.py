from obspy.clients.fdsn import Client
from eqcorrscan.utils.catalog_utils import filter_picks
import numpy as np
import obspy
import matplotlib.pyplot as plt
import eqcorrscan
from obspy.clients.fdsn import Client
from datetime import datetime
import pandas as pd
import dask

start = datetime.now()

# Read in templates
tribe = eqcorrscan.core.match_filter.tribe.read_tribe('growclust_templates_sep2017.tgz')

print('Read in the templates!')

# Let's just use the first template as a test
# tribe = eqcorrscan.core.match_filter.tribe.Tribe(templates = tribe[0])

# Now detect!
client = Client('IRIS')
t1 = datetime(2019,9,1)
t2 = datetime(2019,9,3)
time_bins = pd.to_datetime(np.arange(t1,t2,pd.Timedelta(1,'days')))

@dask.delayed
def loop_days(tribe,day1):
    t1 = obspy.UTCDateTime(day1)
    t2 = obspy.UTCDateTime(day1 + pd.Timedelta(1,'days'))
    return tribe.client_detect(client,t1,t2,threshold=6, threshold_type='MAD',trig_int=1,save_progress=False,process_cores=1,ignore_bad_data=True,concurrency='multiprocess')


lazy_results = [loop_days(tribe,time) for time in time_bins]
    
results = dask.compute(lazy_results)

# Combine the results from separate days into one party
families = []
num_templates = len(results[0][0].families)
for i in range(num_templates):
    fam_list = [p[i] for p in results[0]]
    fam = fam_list[0]
    for f in fam_list[1:]:
        fam = fam.append(f)
    families.append(fam)
party = eqcorrscan.core.match_filter.party.Party(families=families)


party.write('detections_sep2017_dask')

print(datetime.now()-start)