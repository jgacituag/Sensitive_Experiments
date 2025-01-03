import csv
# Read the CSV file
import pickle
import numpy as np
import os 
import wrf_module as wrf
import conf   #Load the default configuration
#from multiprocessing import Pool

def run_wrf(args):
    #############################################################################################
    # CONFIGURATION PARAMETERS (in a future version this can be an input to a function)
    #############################################################################################
    row, i = args
    param1, param2, param3, param4 = map(float, row)
    conf = conf.conf
    conf['expname']   = 'TEST_Multi_4_vars'
    conf['run_num']   = i
    conf['modelpath'] = '/vol0003/hp150019/data/jgacitua/Sensitivity_Experiments/em_quarter_ss'
    conf['datapath']  = '/vol0003/hp150019/data/jgacitua/Sensitivity_Experiments'
    conf['nthreads']  = 3 #Number of threads that will be used for running WRF    

    #Parameters controling the shape of the wind profile.
    conf['modify_wind_profile'] = True  #Are we going to modify the original wind profile?
    conf['remove_mean_wind']    = True  #True for convection experiments, false for mountain wave experiments.
    conf['shear_type']          = 'Curved' #Posible choices: 'Linear','Curved'

    #For the Quarter shear case
    conf['total_shear_depth']  = 4000+(param1*(8000-4000))    #Depth of circular shear. (tipical range [4000 , 8000]) 
    conf['int_total_shear']    = 0+(param2*(30-0))      #Total integrated shear (in m/s, tacking curvature into account) (tipical range [0 , 30])
    conf['curved_shear_per']   = param3       #Which proportion of the total shear depth will be curved. (tipical range [0 , 1])

    #For the Linear shear case
    conf['shear_depth_u']    = 0.0      #Heigth of no shear (for the u component a linear type shear is assumed)
    conf['shear_strength_u'] = 5.0e-3   #Shear in m/s^2 (for the u component)
    conf['shear_depth_v']    = 8000.0   #Heigth of no shear (for the v component)
    conf['shear_strehgth_v'] = 0.0      #Maximum V-wind in m/s (for the v component a sine type shear is assumed)
    conf['shear_freq_v']     = 0.5*np.pi    #This controls the wind gire, 0.5pi means quarter circle pi means half a circle, 2pi means full circle)
    conf['shear_amp_v']      = 0.0      #Maximum V-wind in m/s (for the v component a sine type shear is assumed)
    conf['surf_u'] = 0.0   #u-wind component at the surface (only meaningul is some type of surface drag is active)
    conf['surf_v'] = 0.0   #v-wind component at the surface (only meaningul is some type of surface drag is active)
    conf['llj_amp']         = 0.0      #Low level jet local maximum (for the u component)
    conf['llj_h']           = 1500.0   #Low level jet heigth (for the u component)
    conf['llj_width']       = 500.0    #Low level jet width  (for the u component)
    conf['llj_dir']         = 360.0    #Wind direction for the LLJ. 

    #Parameters controling the stability (temperature profile)
    conf['modify_stability'] = True          #Are we going to modify the original stability?
    conf['stability_factor'] = -1.5+(param4*(1.5--1.5))          #Factor controling the stability change (tipical range [-1.5 , 1.5] )
    conf['stability_factor_height'] = 5000 #2500+(param4*(10000-2500)) #Height of maximum warming / cooling. (tipical range [2500 , 10000 ] )

    #Parameters controling the shape of the potential temperature profile.
    conf['modify_theta_profile'] = False  #Are we going to modify the original temperature?
    conf['surf_theta']           = 288.0 #Potential temperature at the lowest model level.
    #Temperature profile will be computed as costant d(theta)/dz in user defined layers.
    conf['theta_layer_limits']  = np.array([0.0,3000.0,15000.0])  #Layer limits.
    conf['dthetadz_layer']      = np.array([4.0,4.0,4.0])/1000.0 #dthetadz corresponding to each layer.

    #Parameters controling the shape of the moisture profile.
    conf['modify_moisture_profile'] = False       #Are we going to modify the original moisture profile?
    conf['dry_run']                 = False       #Assume 0.0 moisture content at all levels?
    conf['low_level_moisture_height'] = 2000.0    #Low level moisture modification will take effect below this level.
    conf['low_level_moisture_mult_factor'] = 0.0  #Moisture modification for low levels (tipical range [-15.0 , 15.0] )   
    conf['mid_level_moisture_height'] = 2000.0    #Mid level moisture modification will take effect above this level.
    conf['mid_level_moisture_mult_factor'] = 0.0  #Moisture modification factor for mid levels (tipical range [-10.0 , 10.0 ] )

    #Parameters controlling the sensitivity function to be used
    conf['sf_type'] = 'mean_precipitation' 
    conf['sf_percentile'] = 99.0 

    conf['run_model'] = True#Are we going to run the model?
    conf['plot_exp']  = False  #Are we going to do detailed plots of the model solution?

    #Parameters controling the model parameters (parameters controlling the model integration)
    conf['model_xdomain'] = 2000.0 * 42  #Do not change (total extension of the domain in the x direction, meters)
    conf['model_nx'] = 80                #Number of grid points in the X direction (integer) 
    conf['model_dx'] = 2000.0            #Model resolution 
    conf['model_dt'] = 12                #Integration time step (seconds, integer)            
    conf['model_nz'] = 41                #Number of vertical levels (integer)
    conf['model_dy'] = conf['model_dx']  #Model resolution in the y-direction (for the moment we kept this equal to the resolution in the x-direction). 
    conf['model_ny'] = conf['model_nx']  #Number of grid points in the y direction (for the moment we keep this constant). 
    conf['model_dt_fract_num'] = 0       #Numerator for fractional time step (from 1 to 10)
    conf['model_dt_fract_den'] = 1       #Denominator for fractional time step (for the moment we kept this equal to 10)

    #############################################################################################
    #  CALL WRF MODEL
    #############################################################################################
    print('start running wrf')
    wrf.run_wrf( conf )


csv_file_path = 'sampling_batch_test.csv'
i=1
with open(csv_file_path, newline='') as csvfile:
    print(f'reading file {i}')
    csv_reader = csv.reader(csvfile)
    next(csv_reader)  # Skip the header row if it exists
    for row_number,row in enumerate(csv_reader):
        run_wrf([row,row_number])

    #pool = Pool(processes=10)
    #pool.map(run_wrf, csv_reader)
    #pool.map(run_wrf, [(row, row_number) for row_number, row in enumerate(csv_reader)])
    #pool.close()
    #pool.join()

    i+=1

