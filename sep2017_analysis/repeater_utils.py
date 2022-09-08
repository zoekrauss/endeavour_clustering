import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from obspy.core.event import read_events
import obspy
from obspy.io.nlloc.core import read_nlloc_hyp
import os
import eqcorrscan
import eqcorrscan.core.match_filter.template
from obspy.clients.fdsn.client import Client
from eqcorrscan.core.template_gen import template_gen
from eqcorrscan import Tribe
import datetime
import ipympl
import matplotlib
import math
import pymap3d
import itertools
import dask
client = Client('IRIS')

import warnings
warnings.filterwarnings('ignore')

def xy2map(xlc,ylc,depth):
    """
    Converts the cartesian coordinates to lat, lon coordinates following the specific geometry of the ETOMO experiment grids
    Where inp is a list containing [x,y,z] of a NLLoc hypocenter
    """

   
    
    # Need to convert coordinates from kilometers to meters in order to agree with pymap3d functions
    xlc = (xlc-60) * 1000
    ylc = (ylc-90) * 1000

    # Values copied from srGeometry
    rot = -21 
    lon0 = -129.0776
    lat0 = 47.9826


    rad = np.pi/180
    csr = np.cos(rot*rad)
    snr = np.sin(rot*rad)
    x   = xlc*csr - ylc*snr
    y   = xlc*snr + ylc*csr

    alt     = 0
    lat0rad = math.radians(lat0)
    lon0rad = math.radians(lon0)

    ell=pymap3d.ellipsoid.Ellipsoid(model='wgs84')
    [lat,lon,z]=pymap3d.enu2geodetic(x,y,0,lat0rad,lon0rad,alt,ell=ell,deg=False)
    lat = math.degrees(lat)
    lon = math.degrees(lon)

    return lon,lat,depth/1000

def map2xy(lat,lon):
    """
    Converts lat, lon coordinates to cartesian coordinates following the specific geometry of the ETOMO experiment grids
    """
    
    # Values copied from srGeometry
    rot = -21 
    lon0 = -129.0776
    lat0 = 47.9826


    ell=pymap3d.ellipsoid.Ellipsoid(model='wgs84')
    [dx,dy,u] = pymap3d.geodetic2enu(lat, lon, 0, lat0, lon0,0, ell=ell, deg=True)

    rota  = rot*np.pi/180;
    sinrota = np.sin(rota);
    cosrota = np.cos(rota);

    x =  cosrota*dx + sinrota*dy;
    y = -sinrota*dx + cosrota*dy;

    x = x/1000
    y = y/1000
    
    return x,y

def nlloc_loc(catalog):
    """
    Function to locate the picked events in an obspy catalog object using NonLinLoc
    
    Note: Picks in the input catalog must have Phase hints ('P' or 'S')
    
    Outputs a new obspy catalog object with an origin class containing the NLLoc location
    """
    
    
    picks = [ev.picks for ev in catalog]
    
    # Write phase file
    filename = '/Users/zoekrauss/NLLOC/Repeaters/phasefiles/test_phase.txt'
    with open(filename, 'w+b') as f: # Open file in write + append mode in binary
        for i,t in enumerate(catalog):
            catalog[i].write(f,'NLLOC_OBS') # Write one event at a time
            f.write('\n'.encode())
    
    # Perform NLLoc location from shell
    os.system("cd /Users/zoekrauss/NLLOC/Repeaters; ./NLLoc ./run/Repeat_freeZ_control.in")
  
    # Read in results
    fname = '/Users/zoekrauss/NLLOC/Repeaters/freeZ_output/freeZ.sum.grid0.loc.hyp'
    loc = read_nlloc_hyp(fname,coordinate_converter=xy2map)
    
    # Add picks to location structure
    # Note: this relies on both structures being sorted ascending in time
    for i,ev in enumerate(loc):
        ev.resource_id = str(1000 + i)
        ev.picks = catalog[i].picks
        
    return loc

def growclust_all(cat,pathname):
    
    """
    Pathname is string of path to directory containing mseeds of events in the catalog
    """
    
    # Write and save event file
    make_growclust_evfile(cat)
    
    # Create cross correlation file
    make_growclust_ccfile(cat,pathname)
   
    # Perform GrowClust relocation from shell and retrieve outputs
    new_cat = run_growclust()
    
    return new_cat
    
def run_growclust():
    
    # Perform GrowClust relocation from shell
    os.system("cd /Users/zoekrauss/GrowClust3D/endeavour; julia run_growclust_nllgrid3D.jl endeavour.grid3D.inp")
    
    # Read in results
    outfile = '/Users/zoekrauss/GrowClust3D/endeavour/data/out/end_out.grid3D.cat'
    new_cat = pd.read_csv(outfile,header=None,delim_whitespace=True)
    new_cat.rename(columns={0:'yr',1:'mon',2:'day',3:'hr',4:'min',5:'sec',6:'evid',7:'latR',8:'lonR',9:'depR',\
                       10:'mag',11:'ev_serialID',12:'clusterID',13:'clustersize',14:'npair',15:'ndiffP',\
                       16:'ndiffS',17:'rmsP',18:'rmsS',19:'erh',20:'erz',21:'ert',22:'latC',23:'lonC',24:'depC'},inplace=True)

    # Make a datetime column
    new_cat = new_cat.astype({'yr':'int','mon':'int','day':'int','hr':'int','min':'int'})
    new_cat['eventTime'] = new_cat[["yr","mon","day","hr","min"]].apply(lambda x: '-'.join(x.values.astype(int).astype(str)), axis="columns")
    new_cat['eventTime'] = new_cat[['eventTime','sec']].apply(lambda x: '-'.join(x.values.astype(str)), axis="columns")
    new_cat['eventTime']=pd.to_datetime(new_cat['eventTime'],format='%Y-%m-%d-%H-%M-%S.%f')
    
    return(new_cat)

def make_growclust_evfile(cat):
    
     # Write and save event file
    origin_ids = [str(event.resource_id) for event in cat]
    origins = [p.origins[0] for p in cat.events]
    mags = np.zeros(len(cat))
    years = [p.time.year for p in origins]
    months = [p.time.month for p in origins]
    days = [p.time.day for p in origins]
    hours = [p.time.hour for p in origins]
    mins = [p.time.minute for p in origins]
    secs = [p.time.second + p.time.microsecond/1000000 for p in origins]
    lats = [p.latitude for p in origins]
    lons = [p.longitude for p in origins]
    depths = [p.depth for p in origins]
    erh = np.zeros(len(cat))
    erz = np.zeros(len(cat))
    rms = np.zeros(len(cat))
    evtext_dict ={'yr':years,'mon':months,'day':days,'hr':hours,'min':mins,'sec':secs,'lat':lats,'lon':lons,'dep':depths,'mag':mags,'errh':erh,'errz':erz,'rms':rms,'evid':origin_ids}
    events = pd.DataFrame.from_dict(evtext_dict)
    ev_fname = '/Users/zoekrauss/GrowClust3D/endeavour/data/in/evlist.txt'
    events.to_csv(ev_fname,sep='\t',header=False,index=False,float_format="%.3f")
    
    return

def make_growclust_ccfile(cat,pathname):
    """
    Pathname is string of path to directory containing mseeds of events in the catalog
    """
    
    origin_ids = [str(event.resource_id) for event in cat]
    
    # Set up looping, so that we get all possible pairs of events but NO repeats!
    loop_ind = list(itertools.combinations(range(len(cat)),2))
    
    @dask.delayed
    def loop_pairs(pair):
        test = event_corr(pair,cat,origin_ids,pathname)
        return test
    
    lazy_results = [loop_pairs(pair) for pair in loop_ind]
    
    results = dask.compute(lazy_results)
    
    with open('/Users/zoekrauss/GrowClust3D/endeavour/data/in/xcordata.txt','w') as fl:
    
        for arr in results[0]:

            # Write event pair ID line
            fl.write('#'+'\t'+arr[0][0]+'\t'+arr[0][1]+'\t'+'0.000')
            fl.write('\n')

            for pick in arr[1]:
                # Write station- and phase-specific line
                fl.write(pick[0]+'\t'+'{:.4f}'.format(pick[1])+'\t'+'{:.4f}'.format(pick[2])+'\t'+pick[3])
                fl.write('\n')
    
    return
    
def pick_corr(phase,sta,chan,t1,t2,ot1,ot2):
    """
    Cross correlates two phase picks from the same channel and calculates their differential travel times
    This function pulls in the waveform data directly from Iris
    
    INPUT:
    phase: 'P' or 'S'; this determines the size of the window cut around the pick time
    sta: string, station code 
    chan: string, channel code
    t1: pick time of first pick as UTCDatetime
    t2: pick time of second pick as UTCDatetime
    ot1: origin time of the earthquake the first pick is for in UTCDatetime
    ot2: origin time of the earthquake the second pick is for in UTCDatetime
    
    RETURNS:
    value: maximum cross correlation value of the waveform cut around the two picks
    dtt: differential travel times of the two picks
    """
    
    if phase=='P':
        t_off=[-0.3,0.5]
        t_off=[-1,1.5]
    if phase=='S':
        t_off = [-1,2.5]
        
    # Pull in waveforms cut around picks
    tr1 = client.get_waveforms('NV',sta,'',chan,t1+t_off[0],t1+t_off[1])
    tr2 = client.get_waveforms('NV',sta,'',chan,t2+t_off[0],t2+t_off[1])
    tr1.filter('bandpass',freqmin=5,freqmax=15)
    tr2.filter('bandpass',freqmin=5,freqmax=15)
    tr1.plot
    tr2.plot
    
    # Cross correlate waveforms, allowing them to shift relative to each other for the case of faulty picks
    xcor = obspy.signal.cross_correlation.correlate(tr1[0],tr2[0],100)
    shift,value = obspy.signal.cross_correlation.xcorr_max(xcor)
    
    # Differential travel times, calculated as tt1-tt2
    dtt = (t1-ot1) - (t2-ot2)
    
    return(value,dtt)


def event_corr(pair,cat,origin_ids,path):
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
    
    
    i = pair[0]
    j = pair[1]

    ot1 = cat[i].origins[0].time
    ot2 = cat[j].origins[0].time
    sta1 = [p.waveform_id.station_code for p in cat[i].picks]
    sta2 = [p.waveform_id.station_code for p in cat[j].picks]
    # pha1 = [a.phase for a in origins[i].arrivals]
    # pha2 = [a.phase for a in origins[j].arrivals]
    pha1 = [p.phase_hint for p in cat[i].picks]
    pha2 = [p.phase_hint for p in cat[j].picks]
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

def download_mseeds(cat,pathname):
    
    """ 
    Pathname must be a string
    """

    
    # Make sure path exists
    if ~os.path.exists(pathname):
        os.mkdir(pathname)
        
        
    for event in cat:

        evid = str(event.resource_id)

        t1 = min([p.time for p in event.picks])-5
        t2 = max([p.time for p in event.picks])+5

        sta = ",".join(np.unique([p.waveform_id.station_code for p in event.picks]))
        chan = ",".join(np.unique([p.waveform_id.channel_code for p in event.picks]))

        st = client.get_waveforms('NV',sta,'',chan,t1,t2)

        st.write(pathname+str(evid)+'.mseed')
    
    return