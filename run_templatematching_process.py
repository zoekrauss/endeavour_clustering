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
import tarfile
import tempfile
from eqcorrscan.core.match_filter.tribe import Tribe
from eqcorrscan.core.match_filter.family import _write_family


def day_detect(day1,tribe):
    
    t1 = obspy.UTCDateTime(day1)
    t2 = obspy.UTCDateTime(day1 + pd.Timedelta(1,'day'))
    client = Client('IRIS')

    try:
        party=tribe.client_detect(client,t1,t2,threshold=6, threshold_type='MAD',trig_int=1,save_progress=False,ignore_bad_data=True,parallel_process=False)
        fbase = 'detections_sep2017/'
        fname = fbase + 'party_' + t1.strftime('%m_%d_%Y')

        # Make sure file can be overwritten
        if os.path.exists(fname+'.tgz'):
            os.remove(fname+'.tgz')

        # Write party to tar file manually
        # Can't use the built-in function from eqcorrscan when running in multiprocess parallel
        # because the way they do it is by making a temporary directory, 
        # writing the tribes individually, and then zipping them up
        # So instead, I riffed on their source code but made sure that every process names
        # the temporary directory uniquely
        with tempfile.TemporaryDirectory('party_' + t1.strftime('%m_%d_%Y')) as temp_dir:
            Tribe([f.template for f in party.families]).write(
                filename=temp_dir, compress=False,
                catalog_format="QUAKEML")
            for i, family in enumerate(party.families):
                name = family.template.name + '_detections.csv'
                name_to_write = os.path.join(temp_dir, name)
                _write_family(family=family, filename=name_to_write)
            with tarfile.open(fname+'.tgz', "w:gz") as tar:
                tar.add(temp_dir, arcname=os.path.basename(fname+'.tgz'))
    
    except:
        print('Matched filter failed for '+ t1.strftime('%m_%d_%Y'))
    
    
    
    



def main():

    start = timer()

    t1 = datetime(2019,9,1)
    t2 = datetime(2019,9,5)
    time_bins = pd.to_datetime(np.arange(t1,t2,pd.Timedelta(1,'days')))

    tribe = eqcorrscan.core.match_filter.tribe.read_tribe('growclust_templates_sep2017.tgz')
    
    
    for i in range(len(time_bins)):
        p=multiprocessing.Process(target=day_detect,args=(time_bins[i],tribe,))
        # Pin created processes in a round-robin following https://sorami-chi.hateblo.jp/entry/2016/04/29/000000
        os.system("taskset -p -c %d %d" % ((i % os.cpu_count()), os.getpid()))
        p.start()
    
    end = timer()

    print(f'elapsed time: {end - start}')


if __name__ == '__main__':
    main()

    