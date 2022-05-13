# Import modules

import numpy as np
import pandas as pd
import shutil
import os
from zipfile import ZipFile
import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import obspy

# Download the mseeds, which is the actual seismic time series data, for desired time period
from EQTransformer.utils.downloader import downloadMseeds

stime = "2017-01-01"
ftime = "2018-01-01"

downloadMseeds(client_list=["IRIS"], stations_json='end_stationlist.json', output_dir="/home/adminuser/lynx/downloads_mseeds",min_lat =47.5 , max_lat = 48.5,min_lon=-129.4,max_lon=-128.8,start_time=stime,end_time=ftime, chunk_size=1,channel_list=["HH[ZNE]","EH[ZNE]", "HH[Z21]","EH[Z21]", "CH[ZNE]","HH[123]"],n_processor=1)


# Perform detections
from EQTransformer.core.mseed_predictor import mseed_predictor
    
mseed_predictor(input_dir='/home/adminuser/lynx/downloads_mseeds',   
                    input_model='EqT_model.h5',
                    stations_json='end_stationlist.json',
                    output_dir='/home/adminuser/lynx/end_2017',
                    detection_threshold=0.2,                
                    P_threshold=0.1,
                    S_threshold=0.1,  
                    number_of_plots=0,
                    batch_size=20,
                    overlap=0.3,
                    output_probabilities=True)

