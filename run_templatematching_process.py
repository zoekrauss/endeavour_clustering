# This follows python multiprocessing pool method from:
# https://zetcode.com/python/multiprocessing/

# This allows for better management than dask-
# Here I've added in a Keyboard Interrupt exception that will correctly close all processes if needed

from obspy.clients.fdsn import Client
import numpy as np
import obspy
import eqcorrscan
from obspy.clients.fdsn import Client
from datetime import datetime
import time
from timeit import default_timer as timer
import pandas as pd
from multiprocessing import Pool, cpu_count, Process
import os
import multiprocessing

def day_detect(day1):
    tribe = eqcorrscan.core.match_filter.tribe.read_tribe('growclust_templates_sep2017.tgz')
    t1 = obspy.UTCDateTime(day1)
    t2 = obspy.UTCDateTime(day1 + pd.Timedelta(1,'days'))
    client = Client('IRIS')
    try:
        party=tribe.client_detect(client,t1,t2,threshold=6, threshold_type='MAD',trig_int=1,save_progress=False,ignore_bad_data=True,parallel_process=False)
        fbase = 'detections_sep2017/'
        fname = fbase + 'party_' + t1.strftime('%m_%d_%Y')

        # Make sure file can be overwritten
        if os.path.exists(fname+'.tgz'):
            os.remove(fname+'.tgz')


        party.write(fname)
    except:
        print('Matched filter failed for '+ t1.strftime('%m_%d_%Y'))


def main():

    start = timer()

    t1 = datetime(2017,9,1)
    t2 = datetime(2018,9,30)
    time_bins = pd.to_datetime(np.arange(t1,t2,pd.Timedelta(1,'days')))

    for i in range(len(time_bins)):
        p=multiprocessing.Process(target=day_detect,args=(time_bins[i],))
        # Pin created processes in a round-robin following https://sorami-chi.hateblo.jp/entry/2016/04/29/000000
        print(i % os.cpu_count())
        print(os.cpu_count())
        print(os.getpid())
        os.system("taskset -p -c %d %d" % ((i % os.cpu_count()), os.getpid()))
        p.start()
    
    end = timer()

    print(f'elapsed time: {end - start}')


if __name__ == '__main__':
    main()
