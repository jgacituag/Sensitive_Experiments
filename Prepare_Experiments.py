import os
import numpy as np
import pickle

def edit_namelist(file_path, conf):
    with open(file_path + '.tmp', 'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                if '@@DT@@' in line:
                    new_file.write(line.replace('@@DT@@', str(conf['model_dt'])))
                elif '@@DT_FRACT_NUM@@' in line:
                    new_file.write(line.replace('@@DT_FRACT_NUM@@', str(conf['model_dt_fract_num'])))
                elif '@@DT_FRACT_DEN@@' in line:
                    new_file.write(line.replace('@@DT_FRACT_DEN@@', str(conf['model_dt_fract_den'])))
                elif '@@DX@@' in line:
                    new_file.write(line.replace('@@DX@@', str(conf['model_dx'])))
                elif '@@DY@@' in line:
                    new_file.write(line.replace('@@DY@@', str(conf['model_dy'])))
                elif '@@NX@@' in line:
                    new_file.write(line.replace('@@NX@@', str(conf['model_nx'])))
                elif '@@NY@@' in line:
                    new_file.write(line.replace('@@NY@@', str(conf['model_ny'])))
                elif '@@NZ@@' in line:
                    new_file.write(line.replace('@@NZ@@', str(conf['model_nz'])))
                else:
                    new_file.write(line)
    os.system('rm -f ' + file_path)
    os.system('mv ' + file_path + '.tmp ' + file_path)

def read_input_sounding(my_file):

    cp=1004.0
    Rd=287.0
    g=9.8
    P0=100000.0
    input_sounding = dict()
    # Read input_sounding
    with open(my_file) as f:
        lines = f.readlines()
    lines = [line[:-1] for line in lines]  # Remove the new line character
    line_split = lines[0].split()
    input_sounding['surf_pressure'] = float(line_split[0])
    input_sounding['surf_theta'] = float(line_split[1])
    input_sounding['surf_qv'] = float(line_split[2])
    input_sounding['nlevs'] = len(lines) - 1
    input_sounding['height'] = np.zeros(input_sounding['nlevs'])
    input_sounding['theta'] = np.zeros(input_sounding['nlevs'])
    input_sounding['qv'] = np.zeros(input_sounding['nlevs'])
    input_sounding['u'] = np.zeros(input_sounding['nlevs'])
    input_sounding['v'] = np.zeros(input_sounding['nlevs'])
    input_sounding['p'] = np.zeros(input_sounding['nlevs'])
    tmp_press = np.zeros( input_sounding['nlevs'] )
    tmp_press[0] = (input_sounding['surf_pressure']*100) ** ( Rd/cp )

    for ii in range(input_sounding['nlevs']):
        line_split = lines[ii + 1].split()
        input_sounding['height'][ii] = float(line_split[0])
        input_sounding['theta'][ii] = float(line_split[1])
        input_sounding['qv'][ii] = float(line_split[2])
        input_sounding['u'][ii] = float(line_split[3])
        input_sounding['v'][ii] = float(line_split[4])

        if ( ii >= 1 ) :
           theta_mean = 0.5 * ( input_sounding['theta'][ii] + input_sounding['theta'][ii-1] )
           dz         = input_sounding['height'][ii] - input_sounding['height'][ii-1]
           tmp_press[ii] = tmp_press[ii-1] - g * ( P0 ** ( Rd/cp ) ) * ( dz ) / ( cp * theta_mean )
    input_sounding['p'] = tmp_press ** ( ( cp / Rd ) ) / 100.0
    input_sounding['t'] = input_sounding['theta'] * tmp_press / ( P0 ** ( Rd/cp ) ) 

    return input_sounding

def write_input_sounding(my_file, input_sounding):
    with open(my_file, 'w') as f:
        f.write(str(input_sounding['surf_pressure']) + ' ' + str(input_sounding['surf_theta']) + ' ' + str(
            input_sounding['surf_qv']))
        f.write('\n')
        for ii in range(input_sounding['nlevs']):
            f.write(str(input_sounding['height'][ii]) + ' ' + str(input_sounding['theta'][ii]) + ' ' + str(
                input_sounding['qv'][ii]) + ' ' + str(input_sounding['u'][ii]) + ' ' + str(input_sounding['v'][ii]))
            f.write('\n')

def modify_input_sounding(my_file, conf):
    # Read input sounding.
    input_sounding = read_input_sounding(my_file)

    ########################################################################
    # Modify the input sounding wind profile.
    ########################################################################
    if conf['modify_wind_profile'] :
      input_sounding['u'][:] = 0.0 #Overwrite the original wind profile
      input_sounding['v'][:] = 0.0

      if conf['shear_type'] == 'Linear' :
         #Linear shear component - U wind
         max_u = conf['shear_depth_u'] * conf['shear_strength_u']
         input_sounding['u'] = input_sounding['height'] * conf['shear_strength_u']
         input_sounding['u'][ input_sounding['height'] > conf['shear_depth_u'] ] = max_u
         #Linear shear component - V wind
         max_v = conf['shear_depth_v'] * conf['shear_strength_v']
         input_sounding['v'] = input_sounding['height'] * conf['shear_strength_v']
         input_sounding['v'][ input_sounding['height'] > conf['shear_depth_v'] ] = max_v

      if conf['shear_type'] == 'Curved' :
         #Maximum turning in wind direction 
         #is a function of how much shear is explained by the curved part of the hodograph
         #The maximum turning is pi (half a circle), if half the shear is associated with the 
         #curved part then a quarter circle hodograph is obtained.
         curved_shear_theta = conf['curved_shear_per'] * np.pi
         #Curved hodograp up to certain level, then linear shear.

         #Curved hodograp part
         if conf['curved_shear_per'] > 0.0:
            curved_int_shear =  conf['int_total_shear'] * conf['curved_shear_per']
            curved_dep_shear =  conf['total_shear_depth'] * conf['curved_shear_per']
            curved_amp_shear =  curved_int_shear / curved_shear_theta

            mask = input_sounding['height'] <=  curved_dep_shear
            input_sounding['v'][mask] = -1.0 * curved_amp_shear * np.sin( (input_sounding['height'][mask] / curved_dep_shear )* curved_shear_theta )
            input_sounding['u'][mask] = -1.0 * curved_amp_shear * np.cos( (input_sounding['height'][mask] / curved_dep_shear )* curved_shear_theta )
         else :
            curved_amp_shear = 0.0

         #Linear shear part
         u_ini = -1.0 * curved_amp_shear * np.cos( curved_shear_theta )  #U and V wind components at the base of the linear shear layer.
         v_ini = -1.0 * curved_amp_shear * np.sin( curved_shear_theta )

         linear_int_shear = conf['int_total_shear'] * (1.0 - conf['curved_shear_per'])

         u_fin = u_ini + linear_int_shear * np.sin( curved_shear_theta )
         v_fin = v_ini - linear_int_shear * np.cos( curved_shear_theta )

         mask = ( input_sounding['height'] >  curved_dep_shear ) & ( input_sounding['height'] <= conf['total_shear_depth'] )
         input_sounding['v'][mask] = v_ini + ( v_fin - v_ini) * ( input_sounding['height'][mask] - curved_dep_shear ) / ( conf['total_shear_depth'] - curved_dep_shear ) 
         input_sounding['u'][mask] = u_ini + ( u_fin - u_ini) * ( input_sounding['height'][mask] - curved_dep_shear ) / ( conf['total_shear_depth'] - curved_dep_shear )

         #Constant wind part.
         mask = input_sounding['height'] >=  conf['total_shear_depth']
         input_sounding['v'][mask] = v_fin
         input_sounding['u'][mask] = u_fin
         linear_amp_shear = curved_int_shear / curved_shear_theta

      #Low level jet component
      input_sounding['u'] = input_sounding['u'] - conf['llj_amp'] * np.exp( -0.5 * ( conf['llj_h'] - input_sounding['height'] )**2 / conf['llj_width']**2  ) * np.sin( np.deg2rad(conf['llj_dir']) )
      input_sounding['v'] = input_sounding['v'] - conf['llj_amp'] * np.exp( -0.5 * ( conf['llj_h'] - input_sounding['height'] )**2 / conf['llj_width']**2  ) * np.cos( np.deg2rad(conf['llj_dir']) )

      #Add the surface wind speed to the entire profile.
      input_sounding['u'] = input_sounding['u'] + conf['surf_u']
      input_sounding['v'] = input_sounding['v'] + conf['surf_v']

      if conf['remove_mean_wind']:
          # Remove 0-6 km mean wind (this will keep the convection near the center of the
          # domain for a longer time)
          mask = input_sounding['height'] < 6000.0
          mean_u = np.mean(input_sounding['u'][mask])
          mean_v = np.mean(input_sounding['v'][mask])
          input_sounding['u'] = input_sounding['u'] - mean_u
          input_sounding['v'] = input_sounding['v'] - mean_v


    ########################################################################
    # Modify the input sounding theta profile.
    ########################################################################
    if conf['modify_theta_profile']:
        input_sounding['theta'][:] = 0.0
        input_sounding['surf_theta'] = conf['surf_theta']
        input_sounding['theta'][0] = conf['surf_theta']

        index = np.zeros(len(input_sounding['height'])).astype(int)
        for ii in range(len(conf['theta_layer_limits'])):
            mask = input_sounding['height'] >= conf['theta_layer_limits'][ii]
            index[mask] = index[mask] + 1
        index = index - 1
        delta_theta = conf['dthetadz_layer'][index]

        input_sounding['theta'][1:] = conf['surf_theta'] + np.cumsum(
            delta_theta[1:] * (input_sounding['height'][1:] - input_sounding['height'][:-1]))

    ########################################################################
    # Modify stability
    ########################################################################
    if conf['modify_stability']  :
       #factor increases linearly with height. So multiplying factor * temperature can
       #increase or decrease the stability.
       factor = np.sin( input_sounding['height'] * np.pi / (2*conf['stability_factor_height']) )
       factor[ factor < 0.0 ] = 0.0
       factor[ input_sounding['height'] > 2*conf['stability_factor_height'] ] = 0.0
       factor = 1.0 + 0.01 * conf['stability_factor'] * factor
       input_sounding['theta'] = factor * input_sounding['theta']

    ########################################################################
    # Modify the input sounding moisture profile.
    ########################################################################
    if conf['modify_moisture_profile']:
       if conf['dry_run']:
          input_sounding['qv'][:] = 0.0  # Dry simulation.

       #Multiplicative increase/decrease of low level moisture
       tmp_z = (-input_sounding['height'] + conf['low_level_moisture_height'])/250.0
       tmp_factor = 1.0 + 0.01 * conf['low_level_moisture_mult_factor'] / ( 1.0 + np.exp( -tmp_z ) )
       input_sounding['qv'] = input_sounding['qv'] * tmp_factor

       #Multiplicative increase/decrease of mid-upper level moisture.
       tmp_z = ( input_sounding['height'] - conf['mid_level_moisture_height'])/250.0
       tmp_factor = 1.0 + 0.01 * conf['mid_level_moisture_mult_factor'] / ( 1.0 + np.exp( -tmp_z ) )
       input_sounding['qv'] = input_sounding['qv'] * tmp_factor

    ########################################################################
    # Check saturation
    ########################################################################

    epsilon = 0.622
    es = 6.112 * np.exp( 17.67 * ( input_sounding['t'] - 273.16 ) / ( input_sounding['t'] - 273.16 + 243.5 ) )                #[in hPa]
    input_sounding['qvs'] = es * epsilon / ( input_sounding['p'] - ( 1 - epsilon ) * es ) * 1000.0                            #[in g/kg]
    input_sounding['qv'][input_sounding['qv'] > input_sounding['qvs'] ] = input_sounding['qvs'][input_sounding['qv'] > input_sounding['qvs'] ]   #[remove qv values over saturation]


    # Write input_sounding
    write_input_sounding(my_file, input_sounding)

def prepare_wrf(conf):
    print('Preparing WRF Experiment')
    # Create output directory
    expoutpath = conf['datapath'] + '/' + conf['expname'] + '/' + '{:03d}'.format( conf['run_num'] )
    os.makedirs(expoutpath, exist_ok=True)
    # Copy configuration file to the output directory
    with open(expoutpath + '/conf_' + '{:03d}'.format(conf['run_num']) + '.pkl', 'wb') as handle:
        pickle.dump(conf, handle, protocol=pickle.HIGHEST_PROTOCOL)

    tmp_dir = expoutpath + '/tmpdir_' + conf['expname'] + '{:03d}'.format( conf['run_num'] )
    os.makedirs(tmp_dir, exist_ok=True)

    if conf['run_model']:
        # Copy wrf executables and data to the temp directory.
        os.system('cp ' + conf['modelpath'] + '/* ' + tmp_dir + '/')

        #Modify the input sounding 
        input_sounding_file = tmp_dir + '/input_sounding'
        modify_input_sounding(input_sounding_file, conf)

    #############################################################################################
    # STEP 2 - EDIT THE MODEL NAMELIST (MODEL INTEGRATION PARAMETERS AMONG OTHERS)
    #############################################################################################

    if conf['run_model']:
        # Modify the namelist.input
        edit_namelist(tmp_dir + '/namelist.input', conf)

    #############################################################################################
    # STEP 3 - RUN THE MODEL
    #############################################################################################
    outdatafile = expoutpath + '/wrfout_' + '{:03d}'.format(conf['run_num']) + '.nc'

    if conf['run_model']:
        os.chdir(tmp_dir)  # Change working directory to temp_dir.name
        # Run the model (to executables ideal.exe prepares the initial conditions for the model
        # integration, wrf.exe integrates the Navier-Stokes equation from the initial conditions
        # and for the period of the simulation).
        print('Running ideal.exe')
        #os.system('./ideal.exe > out_ideal.log')
        print('Running wrf.exe')
        #os.system('export OMP_NUM_THREADS=' + str(conf['nthreads']) + ';./wrf.exe > out_wrf.log')
        # Move the output (netcdf file containing 4D arrays representing spatio-temporal variation
        # of different physical variables such as temperatur, pressure, wind components, clouds)
        #os.system('mv ./wrfout_d01_0001-01-01_00:00:00 ' + outdatafile)
