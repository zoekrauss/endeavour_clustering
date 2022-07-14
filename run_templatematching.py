from obspy.clients.fdsn import Client
from eqcorrscan.utils.catalog_utils import filter_picks
import numpy as np
import obspy
import matplotlib.pyplot as plt
import eqcorrscan
from obspy.clients.fdsn import Client
from datetime import datetime

start = datetime.now()

# Read in templates
tribe = eqcorrscan.core.match_filter.tribe.read_tribe('growclust_templates_sep2017.tgz')

print('Read in the templates!')

# Let's just use the first template as a test
# tribe = eqcorrscan.core.match_filter.tribe.Tribe(templates = tribe[0])

# Now detect!
client = Client('IRIS')
t1 = obspy.UTCDateTime("2017-09-01")
t2 = obspy.UTCDateTime("2017-09-30")

party = tribe.client_detect(client,t1,t2,threshold=6, threshold_type='MAD',trig_int=1,save_progress=False,process_cores=1,ignore_bad_data=True,save_progress=True)

party.write('detections_sep2017')

print(datetime.now()-start)