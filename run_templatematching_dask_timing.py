# Written to perform template matching using eqcorrscan over the course of a month

# This script is meant to test how parallelizing over different variations in time affects performance

# Should write a text file containing the time it took to do each of the splits specified in time_splits


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

# Read in templates
tribe = eqcorrscan.core.match_filter.tribe.read_tribe('growclust_templates_sep2017.tgz')

print('Read in the templates!')


# Specify time range
client = Client('IRIS')
t1 = datetime(2019,9,1)
t2 = datetime(2019,9,30)


# Specify the time splits you want to test
time_splits = [1,2,3,5,6,10,30] # number of days to split the month up into


# Set up the dask parallelization framework which splits the process over days
@dask.delayed
def loop_days(tribe,day1,split):
    t1 = obspy.UTCDateTime(day1)
    t2 = obspy.UTCDateTime(day1 + pd.Timedelta(split,'days'))
    try:
        test = tribe.client_detect(client,t1,t2,threshold=6,threshold_type='MAD',trig_int=1,save_progress=False,process_cores=1,ignore_bad_data=True,concurrency='multiprocess')
    except:
        print('matched filter failed')
    return t1

# Now loop through the different splits and get the timing for each
durations = []
for split in time_splits:
    print('Starting to loop by ' + str(split) + ' days!!!')
    time_bins = pd.to_datetime(np.arange(t1,t2,pd.Timedelta(split,'days')))
    lazy_results = [loop_days(tribe,time,split) for time in time_bins]
    start = datetime.now()
    results = dask.compute(lazy_results)
    dur = datetime.now()-start
    durations.append(dur)


# Write the results to file
textfile = open("dask_durations.txt", "w")
for element in durations:
    textfile.write(str(element.seconds) + "\n")
textfile.close()