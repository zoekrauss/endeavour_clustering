import numpy as np
import pandas as pd
import shutil
import os
from zipfile import ZipFile
import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import obspy
import glob2 as glob
from obspy.clients.fdsn import Client
from collections import defaultdict
import pyarrow
import itertools
import obspy.signal.cross_correlation
import dask
client = Client('IRIS')

# Get original events
file = 'endquakes_2017.xml'
cat = obspy.core.event.read_events(file)

# FILTER IF SO DESIRED
cat = cat.filter("time > 2017-09-01","time < 2017-09-30","longitude > -129.15","longitude < -129.05","latitude < 48.05","latitude > 47.9")

origins = [p.origins[0] for p in cat.events]
mags = np.empty(len(cat))
for i,ev in enumerate(cat.events):
    if len(ev.magnitudes)>0:
        mags[i] = float(ev.magnitudes[0].mag)
years = [p.time.year for p in origins]
months = [p.time.month for p in origins]
days = [p.time.day for p in origins]
hours = [p.time.hour for p in origins]
mins = [p.time.minute for p in origins]
secs = [p.time.second + p.time.microsecond/1000000 for p in origins]
lats = [p.latitude for p in origins]
lons = [p.longitude for p in origins]
depths = [p.depth for p in origins]
arrivals = [len(p.arrivals) for p in origins]
origin_ids = [str(p.resource_id)[-6:] for p in cat]
erh = [float(p.origin_uncertainty.horizontal_uncertainty) for p in origins]
erz = [float(p.depth_errors['uncertainty']) for p in origins]
rms = np.zeros(len(cat))

# Set up looping, so that we get all possible pairs of events but NO repeats!
loop_ind = list(itertools.combinations(range(len(cat)),2))

def event_corr(pair):
    """
    INPUT:
    A tuple of indices of events in the catalog
    
    For the two events corresponding to those indices, finds the common picks between them and calculates
    their cross correlation coefficient and differential travel times
    
    RETURNS:
    A numpy array designed to be easily written to a text file in the format needed for GrowClust
    Array has format ((ev1id,ev2id),(commonpick1,commonpick2,...))
    Where each commonpick array has the format (station,differential travel time,cross correlation coefficient,phase)
    """
    
    path = 'sep2017_mseed/'
    
    i = pair[0]
    j = pair[1]

    ot1 = cat[i].origins[0].time
    ot2 = cat[j].origins[0].time
    sta1 = [p.waveform_id.station_code for p in cat[i].picks]
    sta2 = [p.waveform_id.station_code for p in cat[j].picks]
    pha1 = [a.phase for a in origins[i].arrivals]
    pha2 = [a.phase for a in origins[j].arrivals]
    picks1 = [a+'_'+pha1[k] for k,a in enumerate(sta1)]
    picks2 = [a+'_'+pha2[k] for k,a in enumerate(sta2)]
    common_picks = list(set(picks1).intersection(set(picks2)))

    pick_arr = np.empty((len(common_picks),4),dtype=object)
    for k,p in enumerate(common_picks):
        # Find corresponding pick info within each event
        # Send to pick_corr to get information
        sta,phase = p.split('_')

        sta1_ind = [ ind for ind in range(len(sta1)) if sta1[ind] == sta ]
        pha1_ind = [ ind for ind in range(len(pha1)) if pha1[ind] == phase ]
        pick1_ind = list(set(sta1_ind).intersection(set(pha1_ind)))
        pick1 = cat[i].picks[pick1_ind[0]]

        sta2_ind = [ ind for ind in range(len(sta2)) if sta2[ind] == sta ]
        pha2_ind = [ ind for ind in range(len(pha2)) if pha2[ind] == phase ]
        pick2_ind = list(set(sta2_ind).intersection(set(pha2_ind)))
        pick2 = cat[j].picks[pick2_ind[0]]

        chan = pick1.waveform_id.channel_code
        xcor_val,dtt = pick_corr_mseed(phase,sta,chan,origin_ids[i],origin_ids[j],pick1.time,pick2.time,ot1,ot2,path)
        # xcor_val,dtt = pick_corr(phase,sta,chan,pick1.time,pick2.time,ot1,ot2)

        # Write station- and phase-specific line
        
        pick_arr[k][0]=sta
        pick_arr[k][1]=dtt
        pick_arr[k][2]=xcor_val
        pick_arr[k][3] = phase
        
    
    return(((origin_ids[i],origin_ids[j]),pick_arr))

def pick_corr_mseed(phase,sta,chan,evid1,evid2,t1,t2,ot1,ot2,path):
    """
    Cross correlates two phase picks from the same channel and calculates their differential travel times
    This function reads waveform data from local miniseed files
    
    INPUT:
    phase: 'P' or 'S'; this determines the size of the window cut around the pick time
    sta: string, station code 
    chan: string, channel code
    t1: pick time of first pick as UTCDatetime
    t2: pick time of second pick as UTCDatetime
    ot1: origin time of the earthquake the first pick is for in UTCDatetime
    ot2: origin time of the earthquake the second pick is for in UTCDatetime
    path: path to the folder that stores the mseed files, as a string. eg 'sep2017_mseed/'
    
    RETURNS:
    value: maximum cross correlation value of the waveform cut around the two picks
    dtt: differential travel times of the two picks
    """
    
    if phase=='P':
        t_off=[-0.2,0.6]
        # t_off=[-1,1.5]
    if phase=='S':
        t_off = [-0.3,0.6]
        
    # Pull in waveforms
    fid1 = path+str(evid1)+'.mseed'
    fid2 = path+str(evid2)+'.mseed'
    tr1 = obspy.read(fid1).select(station=sta,channel=chan).detrend()
    tr2 = obspy.read(fid2).select(station=sta,channel=chan).detrend()
    
    # Filter
    tr1.filter('bandpass',freqmin=8,freqmax=35)
    tr2.filter('bandpass',freqmin=8,freqmax=35)
    
    # Cut the waveforms around picks
    tr1.trim(starttime=t1+t_off[0],endtime=t1+t_off[1])
    tr2.trim(starttime=t2+t_off[0],endtime=t2+t_off[1])
    
    
    # Cross correlate waveforms, allowing them to shift relative to each other for the case of faulty picks
    xcor = obspy.signal.cross_correlation.correlate(tr1[0],tr2[0],100)
    shift,value = obspy.signal.cross_correlation.xcorr_max(xcor)
    
    
    
    # Differential travel times, calculated using the shift
    new_t2 = t2-(shift/tr2[0].stats.sampling_rate)
    dtt = (t1-ot1) - (new_t2-ot2)
    
    return(value,dtt)

@dask.delayed
def loop_pairs(pair):
    try:
        test = event_corr(pair)
        return test
    except: 
        print(pair)


lazy_results = [loop_pairs(pair) for pair in loop_ind]


results = dask.compute(lazy_results)

with open('end_xcordata_0520_2.txt','w') as fl:
    
    for arr in results[0]:
        
        # Write event pair ID line
        fl.write('#'+'\t'+arr[0][0]+'\t'+arr[0][1]+'\t'+'0.000')
        fl.write('\n')
        
        for pick in arr[1]:
            # Write station- and phase-specific line
            fl.write(pick[0]+'\t'+'{:.4f}'.format(pick[1])+'\t'+'{:.4f}'.format(pick[2])+'\t'+pick[3])
            fl.write('\n')