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
tribe = eqcorrscan.core.match_filter.tribe.read_tribe('filtered_templates_july2019.tgz')

# Let's just use the first template as a test
tribe = eqcorrscan.core.match_filter.tribe.Tribe(templates = tribe[0])

# Now detect!
client = Client('IRIS')
t1 = obspy.UTCDateTime("2019-07-01-17")
t2 = obspy.UTCDateTime("2019-07-01-18")

party = tribe.client_detect(client,t1,t2,threshold=8, threshold_type='MAD',trig_int=1)

party.write('test_detections')

print(datetime.now()-start)