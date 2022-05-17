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
                    output_probabilities=True,gpuid=0)

