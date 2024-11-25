from netCDF4 import Dataset
from wrf import (getvar, to_np, interplevel, extract_times)
import glob
import numpy as np
import pickle
import os
import datetime as dt

# Constantes utiles
cp = 1004.0  # Calor especifico del aire seco a presion constante en J/K
Rd = 287.0  # Constante del aire seco
gg = 9.81  # Gravedad en m/s
Lv = 2264.705e3  # Calor latente de vaporizacion en J/Kg

def get_data_vslice(data_path, slice_width=1, slice_index=0, slice_z=None, zres=None, slice_type='vy', t_start=0,
                    t_end=0, force=False):
    if slice_type == 'vy' or slice_type == 'vx':
        file_name = data_path + '/wrfout_' + slice_type + '_slice_index_' + str(slice_index) + '_t_start_' + str(
            t_start) + '_t_end_' + str(t_end) + '.pkl'
    if slice_type == 'h':
        file_name = data_path + '/wrfout_' + slice_type + '_slice_z_' + str(slice_z) + '_t_start_' + str(
            t_start) + '_t_end_' + str(t_end) + '.pkl'

    if slice_type == 'h' and (slice_z is None or zres is None):
        print('Error: Si se hace un corte horizontal hay que especificar slice_z y zres')
        return

    # Vamos a guardar una cross section a x o y constantes. Si fijamos un x o un y, tomamos width puntos alrededor de ese x o y para poder calcular
    # derivadas espaciales centradas en el corte. La idea es tener todos los datos necesarios para poder hacer calculos sobre ese corte en particular (derivadas temporales y espaciales).

    if os.path.exists(file_name) & (not force):
        with open(file_name, 'rb') as handle:
            my_data = pickle.load(handle)

    else:

        my_data = dict()
        my_data['file_name'] = file_name
        my_data['times'] = []

        my_data['file_list'] = glob.glob(data_path + '/wrfout_d01_*')  # Busco todos los wrfout en la carpeta indicada.
        my_data['file_list'].sort()

        # Leo el primer archivo de la lista de donde voy a tomar el perfil de referencia.
        ncfile = Dataset(my_data['file_list'][0])
        tmp = to_np(getvar(ncfile, "P"))
        nx = tmp.shape[2]  # Order is Z=0 , Y=1 , X=0
        ny = tmp.shape[1]
        nz = tmp.shape[0]
        nt = len(my_data['file_list'])

        if (t_start != 0 or t_end != 0) and (t_end > t_start):
            ts = t_start
            te = t_end
        else:
            ts = 0
            te = nt

        if slice_type == 'vy' or slice_type == 'vx':  # Vertical slice
            if slice_type == 'vy':
                ys = slice_index - slice_width
                ye = slice_index + slice_width + 1
                xs = 0
                xe = nx
            if slice_type == 'vx':
                xs = slice_index - slice_width
                xe = slice_index + slice_width + 1
                ys = 0
                ye = ny

            zs = 0
            ze = nz
            snz = ze - zs

        if slice_type == 'h':  # Horizontal slice
            ys = 0
            ye = ny
            xs = 0
            xe = nx
            snz = 2 * slice_width + 1

        snx = xe - xs
        sny = ye - ys
        snt = te - ts

        my_data['nx'] = snx
        my_data['ny'] = sny
        my_data['nz'] = snz
        my_data['nt'] = snt

        # Allocate memory
        my_data['p'] = np.zeros((snz, sny, snx, snt))
        my_data['u'] = np.zeros((snz, sny, snx, snt))
        my_data['v'] = np.zeros((snz, sny, snx, snt))
        my_data['w'] = np.zeros((snz, sny, snx, snt))
        my_data['ref'] = np.zeros((snz, sny, snx, snt))
        my_data['t'] = np.zeros((snz, sny, snx, snt))
        my_data['qv'] = np.zeros((snz, sny, snx, snt))
        my_data['qc'] = np.zeros((snz, sny, snx, snt))
        my_data['qr'] = np.zeros((snz, sny, snx, snt))
        my_data['qi'] = np.zeros((snz, sny, snx, snt))
        my_data['qs'] = np.zeros((snz, sny, snx, snt))
        my_data['qg'] = np.zeros((snz, sny, snx, snt))
        my_data['z'] = np.zeros((snz, sny, snx, snt))
        my_data['thetae'] = np.zeros((snz, sny, snx, snt))
        my_data['theta'] = np.zeros((snz, sny, snx, snt))
        my_data['rh'] = np.zeros((snz, sny, snx, snt))
        my_data['h_diabatic'] = np.zeros((snz, sny, snx, snt))
        my_data['qv_diabatic'] = np.zeros((snz, sny, snx, snt))

        my_data['file_list'] = my_data['file_list'][ts:te]
        # Leo los datos y los guardo en el diccionario del corte.
        for itime in range(ts, te):
            ctime = itime - ts
            print('Reading file ' + my_data['file_list'][ctime])
            ncfile = Dataset(my_data['file_list'][ctime])
            if slice_type == 'vy' or slice_type == 'vx':  # Vertical slice
                my_data['p'][:, :, :, ctime] = to_np(getvar(ncfile, "pres", units='Pa'))[zs:ze, ys:ye, xs:xe]
                my_data['u'][:, :, :, ctime] = to_np(getvar(ncfile, "ua", units='ms-1'))[zs:ze, ys:ye, xs:xe]
                my_data['v'][:, :, :, ctime] = to_np(getvar(ncfile, "va", units='ms-1'))[zs:ze, ys:ye, xs:xe]
                my_data['w'][:, :, :, ctime] = to_np(getvar(ncfile, "wa", units='ms-1'))[zs:ze, ys:ye, xs:xe]
                my_data['ref'][:, :, :, ctime] = to_np(getvar(ncfile, "dbz"))[zs:ze, ys:ye, xs:xe]
                my_data['t'][:, :, :, ctime] = to_np(getvar(ncfile, "temp", units='K'))[zs:ze, ys:ye, xs:xe]
                my_data['qv'][:, :, :, ctime] = to_np(getvar(ncfile, "QVAPOR"))[zs:ze, ys:ye, xs:xe]
                try:
                    my_data['qc'][:, :, :, ctime] = to_np(getvar(ncfile, "QCLOUD"))[zs:ze, ys:ye, xs:xe]
                    my_data['qr'][:, :, :, ctime] = to_np(getvar(ncfile, "QRAIN"))[zs:ze, ys:ye, xs:xe]
                    my_data['qi'][:, :, :, ctime] = to_np(getvar(ncfile, "QICE"))[zs:ze, ys:ye, xs:xe]
                    my_data['qs'][:, :, :, ctime] = to_np(getvar(ncfile, "QSNOW"))[zs:ze, ys:ye, xs:xe]
                    my_data['qg'][:, :, :, ctime] = to_np(getvar(ncfile, "QGRAUP"))[zs:ze, ys:ye, xs:xe]
                except:
                    print('Warning no diabatic data found')

                my_data['z'][:, :, :, ctime] = to_np(getvar(ncfile, "height", units='m'))[zs:ze, ys:ye, xs:xe]
                my_data['thetae'][:, :, :, ctime] = to_np(getvar(ncfile, "theta_e", units='K'))[zs:ze, ys:ye, xs:xe]
                my_data['theta'][:, :, :, ctime] = to_np(getvar(ncfile, "theta", units='K'))[zs:ze, ys:ye, xs:xe]
                my_data['rh'][:, :, :, ctime] = to_np(getvar(ncfile, "rh"))[zs:ze, ys:ye, xs:xe]
                try:
                    my_data['h_diabatic'][:, :, :, ctime] = to_np(getvar(ncfile, "H_DIABATIC"))[zs:ze, ys:ye, xs:xe]
                    my_data['qv_diabatic'][:, :, :, ctime] = to_np(getvar(ncfile, "QV_DIABATIC"))[zs:ze, ys:ye, xs:xe]
                except:
                    print('Warning no diabatic data found')
                my_data['times'].append(
                    dt.datetime.strptime(os.path.basename(my_data['file_list'][ctime])[11:], '%Y-%m-%d_%H:%M:%S'))
            if slice_type == 'h':  # Horizontal slice at constant height
                z_levels = np.arange(slice_z - slice_width * zres, slice_z + 2.0 * slice_width * zres, zres)

                z_ = getvar(ncfile, "height", units='m')
                p_ = getvar(ncfile, "pres", units='Pa')
                u_ = getvar(ncfile, "ua", units='ms-1')
                v_ = getvar(ncfile, "va", units='ms-1')
                w_ = getvar(ncfile, "wa", units='ms-1')
                ref_ = getvar(ncfile, "dbz")
                t_ = getvar(ncfile, "temp", units='K')
                qv_ = getvar(ncfile, "QVAPOR")
                try:
                    qc_ = getvar(ncfile, "QCLOUD")
                    qr_ = getvar(ncfile, "QRAIN")
                    qi_ = getvar(ncfile, "QICE")
                    qs_ = getvar(ncfile, "QSNOW")
                    qg_ = getvar(ncfile, "QGRAUP")
                except:
                    print('Warning no condensate data found')

                thetae_ = getvar(ncfile, "theta_e", units='K')
                theta_ = getvar(ncfile, "theta", units='K')
                rh_ = getvar(ncfile, "rh")
                try:
                    h_diabatic_ = getvar(ncfile, "H_DIABATIC")
                    qv_diabatic_ = getvar(ncfile, "QV_DIABATIC")
                    my_data['h_diabatic'][:, :, :, ctime] = interplevel(h_diabatic_, z_, z_levels)
                    my_data['qv_diabatic'][:, :, :, ctime] = interplevel(qv_diabatic_, z_, z_levels)
                except:
                    print('Warning no diabatic data found ')
                my_data['p'][:, :, :, ctime] = interplevel(p_, z_, z_levels)
                my_data['u'][:, :, :, ctime] = interplevel(u_, z_, z_levels)
                my_data['v'][:, :, :, ctime] = interplevel(v_, z_, z_levels)
                my_data['w'][:, :, :, ctime] = interplevel(w_, z_, z_levels)
                my_data['ref'][:, :, :, ctime] = interplevel(ref_, z_, z_levels)
                my_data['t'][:, :, :, ctime] = interplevel(t_, z_, z_levels)
                my_data['qv'][:, :, :, ctime] = interplevel(qv_, z_, z_levels)
                try:
                    my_data['qc'][:, :, :, ctime] = interplevel(qc_, z_, z_levels)
                    my_data['qr'][:, :, :, ctime] = interplevel(qr_, z_, z_levels)
                    my_data['qi'][:, :, :, ctime] = interplevel(qi_, z_, z_levels)
                    my_data['qs'][:, :, :, ctime] = interplevel(qs_, z_, z_levels)
                    my_data['qg'][:, :, :, ctime] = interplevel(qg_, z_, z_levels)
                except:
                    print('Warning no condensate data found')

                my_data['thetae'][:, :, :, ctime] = interplevel(thetae_, z_, z_levels)
                my_data['theta'][:, :, :, ctime] = interplevel(theta_, z_, z_levels)
                my_data['rh'][:, :, :, ctime] = interplevel(rh_, z_, z_levels)

                my_data['times'].append(
                    dt.datetime.strptime(os.path.basename(my_data['file_list'][ctime])[11:], '%Y-%m-%d_%H:%M:%S'))

        if slice_type == 'h':
            for iz, my_z in enumerate(z_levels):
                my_data['z'][iz, :, :, :] = my_z

        if nt == 1:
            my_data['dt'] = 1.0
        else:
            my_data['dt'] = (my_data['times'][1] - my_data['times'][0]).seconds
        my_data['dx'] = ncfile.DX
        my_data['dy'] = ncfile.DY
        my_data['file_name'] = file_name
        my_data['slice_type'] = slice_type
        my_data['slice_index'] = slice_index
        my_data['slice_width'] = slice_width
        my_data['slice_z'] = slice_z
        my_data['zres'] = zres
        my_data['ts'] = ts
        my_data['te'] = te

        my_data['p'] = np.float32(my_data['p'])
        my_data['u'] = np.float32(my_data['u'])
        my_data['v'] = np.float32(my_data['v'])
        my_data['w'] = np.float32(my_data['w'])
        my_data['ref'] = np.float32(my_data['ref'])
        my_data['t'] = np.float32(my_data['t'])
        my_data['qv'] = np.float32(my_data['qv'])
        my_data['qc'] = np.float32(my_data['qc'])
        my_data['qr'] = np.float32(my_data['qr'])
        my_data['qi'] = np.float32(my_data['qi'])
        my_data['qs'] = np.float32(my_data['qs'])
        my_data['qg'] = np.float32(my_data['qg'])
        my_data['z'] = np.float32(my_data['z'])
        my_data['thetae'] = np.float32(my_data['thetae'])
        my_data['theta'] = np.float32(my_data['theta'])
        my_data['rh'] = np.float32(my_data['rh'])
        my_data['h_diabatic'] = np.float32(my_data['h_diabatic'])
        my_data['qv_diabatic'] = np.float32(my_data['qv_diabatic'])

        with open(file_name, 'wb') as handle:
            pickle.dump(my_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return my_data


def get_data_vslice_single_file(my_file, slice_width=1, slice_index=0, slice_z=None, zres=None, slice_type='vy',
                                t_start=0, t_end=0, force=False, write_pkl=True):
    data_path = os.path.dirname(my_file)
    if slice_type == 'vy' or slice_type == 'vx':
        file_name = data_path + '/wrfout_' + slice_type + '_slice_index_' + str(slice_index) + '_t_start_' + str(
            t_start) + '_t_end_' + str(t_end) + '.pkl'
    if slice_type == 'h':
        file_name = data_path + '/wrfout_' + slice_type + '_slice_z_' + str(slice_z) + '_t_start_' + str(
            t_start) + '_t_end_' + str(t_end) + '.pkl'

    if slice_type == 'h' and (slice_z is None or zres is None):
        print('Error: Si se hace un corte horizontal hay que especificar slice_z y zres')
        return

    # Vamos a guardar una cross section a x o y constantes. Si fijamos un x o un y, tomamos width puntos alrededor de ese x o y para poder calcular
    # derivadas espaciales centradas en el corte. La idea es tener todos los datos necesarios para poder hacer calculos sobre ese corte en particular (derivadas temporales y espaciales).

    if os.path.exists(file_name) & (not force):
        with open(file_name, 'rb') as handle:
            my_data = pickle.load(handle)

    else:

        my_data = dict()
        my_data['file_name'] = file_name
        my_data['times'] = []

        # my_data['file_list'] = glob.glob(data_path + '/wrfout_d01_*')   #Busco todos los wrfout en la carpeta indicada.
        # my_data['file_list'].sort()

        # Leo el primer archivo de la lista de donde voy a tomar el perfil de referencia.
        ncfile = Dataset(my_file)
        tmp = to_np(getvar(ncfile, "P"))
        nx = tmp.shape[2]  # Order is Z=0 , Y=1 , X=0
        ny = tmp.shape[1]
        nz = tmp.shape[0]
        raw_times = extract_times(ncfile, None, meta=False)
        # Convert numpy.datetime64 into datetime format.
        file_times = [dt.datetime.strptime(np.datetime_as_string(tmpdate, timezone='UTC')[0:19], '%Y-%m-%dT%H:%M:%S')
                      for tmpdate in raw_times]
        nt = len(file_times)

        # nt=len( my_data['file_list'] )

        if (t_start != 0 or t_end != 0) and (t_end > t_start):
            ts = t_start
            te = t_end
        else:
            ts = 0
            te = nt

        if slice_type == 'vy' or slice_type == 'vx':  # Vertical slice
            if slice_type == 'vy':
                ys = max(slice_index - slice_width, 0)
                ye = min(slice_index + slice_width + 1, ny)
                xs = 0
                xe = nx
            if slice_type == 'vx':
                xs = max(slice_index - slice_width, 0)
                xe = min(slice_index + slice_width + 1, nx)
                ys = 0
                ye = ny

            zs = 0
            ze = nz
            snz = ze - zs

        if slice_type == 'h':  # Horizontal slice
            ys = 0
            ye = ny
            xs = 0
            xe = nx
            snz = 2 * slice_width + 1

        snx = xe - xs
        sny = ye - ys
        snt = te - ts

        my_data['nx'] = snx
        my_data['ny'] = sny
        my_data['nz'] = snz
        my_data['nt'] = snt

        # Allocate memory
        my_data['p'] = np.zeros((snz, sny, snx, snt))
        my_data['u'] = np.zeros((snz, sny, snx, snt))
        my_data['v'] = np.zeros((snz, sny, snx, snt))
        my_data['w'] = np.zeros((snz, sny, snx, snt))
        my_data['ref'] = np.zeros((snz, sny, snx, snt))
        my_data['t'] = np.zeros((snz, sny, snx, snt))
        my_data['qv'] = np.zeros((snz, sny, snx, snt))
        my_data['qc'] = np.zeros((snz, sny, snx, snt))
        my_data['qr'] = np.zeros((snz, sny, snx, snt))
        my_data['qi'] = np.zeros((snz, sny, snx, snt))
        my_data['qs'] = np.zeros((snz, sny, snx, snt))
        my_data['qg'] = np.zeros((snz, sny, snx, snt))
        my_data['z'] = np.zeros((snz, sny, snx, snt))
        my_data['thetae'] = np.zeros((snz, sny, snx, snt))
        my_data['theta'] = np.zeros((snz, sny, snx, snt))
        my_data['rh'] = np.zeros((snz, sny, snx, snt))
        my_data['h_diabatic'] = np.zeros((snz, sny, snx, snt))
        my_data['qv_diabatic'] = np.zeros((snz, sny, snx, snt))

        # my_data['file_list'] = my_data['file_list'][ts:te]
        # Leo los datos y los guardo en el diccionario del corte.
        for itime in range(ts, te):
            ctime = itime - ts
            print('Reading file ' + my_file)
            ncfile = Dataset(my_file)
            if slice_type == 'vy' or slice_type == 'vx':  # Vertical slice
                my_data['p'][:, :, :, ctime] = to_np(getvar(ncfile, "pres", units='Pa', timeidx=itime))[zs:ze, ys:ye,
                                               xs:xe]
                my_data['u'][:, :, :, ctime] = to_np(getvar(ncfile, "ua", units='ms-1', timeidx=itime))[zs:ze, ys:ye,
                                               xs:xe]
                my_data['v'][:, :, :, ctime] = to_np(getvar(ncfile, "va", units='ms-1', timeidx=itime))[zs:ze, ys:ye,
                                               xs:xe]
                my_data['w'][:, :, :, ctime] = to_np(getvar(ncfile, "wa", units='ms-1', timeidx=itime))[zs:ze, ys:ye,
                                               xs:xe]
                my_data['t'][:, :, :, ctime] = to_np(getvar(ncfile, "temp", units='K', timeidx=itime))[zs:ze, ys:ye,
                                               xs:xe]
                my_data['qv'][:, :, :, ctime] = to_np(getvar(ncfile, "QVAPOR", timeidx=itime))[zs:ze, ys:ye, xs:xe]
                try:
                    my_data['ref'][:, :, :, ctime] = to_np(getvar(ncfile, "dbz", timeidx=itime))[zs:ze, ys:ye, xs:xe]
                    my_data['qc'][:, :, :, ctime] = to_np(getvar(ncfile, "QCLOUD", timeidx=itime))[zs:ze, ys:ye, xs:xe]
                    my_data['qr'][:, :, :, ctime] = to_np(getvar(ncfile, "QRAIN", timeidx=itime))[zs:ze, ys:ye, xs:xe]
                    my_data['qi'][:, :, :, ctime] = to_np(getvar(ncfile, "QICE", timeidx=itime))[zs:ze, ys:ye, xs:xe]
                    my_data['qs'][:, :, :, ctime] = to_np(getvar(ncfile, "QSNOW", timeidx=itime))[zs:ze, ys:ye, xs:xe]
                    my_data['qg'][:, :, :, ctime] = to_np(getvar(ncfile, "QGRAUP", timeidx=itime))[zs:ze, ys:ye, xs:xe]
                except:
                    print('Warning no condensate data found')
                my_data['z'][:, :, :, ctime] = to_np(getvar(ncfile, "height", units='m', timeidx=itime))[zs:ze, ys:ye,
                                               xs:xe]
                my_data['thetae'][:, :, :, ctime] = to_np(getvar(ncfile, "theta_e", units='K', timeidx=itime))[zs:ze,
                                                    ys:ye, xs:xe]
                my_data['theta'][:, :, :, ctime] = to_np(getvar(ncfile, "theta", units='K', timeidx=itime))[zs:ze,
                                                   ys:ye, xs:xe]
                my_data['rh'][:, :, :, ctime] = to_np(getvar(ncfile, "rh", timeidx=itime))[zs:ze, ys:ye, xs:xe]
                try:
                    my_data['h_diabatic'][:, :, :, ctime] = to_np(getvar(ncfile, "H_DIABATIC", timeidx=itime))[zs:ze,
                                                            ys:ye, xs:xe]
                    my_data['qv_diabatic'][:, :, :, ctime] = to_np(getvar(ncfile, "QV_DIABATIC", timeidx=itime))[zs:ze,
                                                             ys:ye, xs:xe]
                except:
                    print('Warning no diabatic data found')
                my_data['times'].append(file_times[itime])
            if slice_type == 'h':  # Horizontal slice at constant height
                z_levels = np.arange(slice_z - slice_width * zres, slice_z + 2.0 * slice_width * zres, zres)

                z_ = getvar(ncfile, "height", units='m', timeidx=itime)
                p_ = getvar(ncfile, "pres", units='Pa', timeidx=itime)
                u_ = getvar(ncfile, "ua", units='ms-1', timeidx=itime)
                v_ = getvar(ncfile, "va", units='ms-1', timeidx=itime)
                w_ = getvar(ncfile, "wa", units='ms-1', timeidx=itime)
                t_ = getvar(ncfile, "temp", units='K', timeidx=itime)
                qv_ = getvar(ncfile, "QVAPOR", timeidx=itime)
                try:
                    ref_ = getvar(ncfile, "dbz", timeidx=itime)
                    qc_ = getvar(ncfile, "QCLOUD", timeidx=itime)
                    qr_ = getvar(ncfile, "QRAIN", timeidx=itime)
                    qi_ = getvar(ncfile, "QICE", timeidx=itime)
                    qs_ = getvar(ncfile, "QSNOW", timeidx=itime)
                    qg_ = getvar(ncfile, "QGRAUP", timeidx=itime)
                    my_data['ref'][:, :, :, ctime] = interplevel(ref_, z_, z_levels)
                    my_data['qc'][:, :, :, ctime] = interplevel(qc_, z_, z_levels)
                    my_data['qr'][:, :, :, ctime] = interplevel(qr_, z_, z_levels)
                    my_data['qi'][:, :, :, ctime] = interplevel(qi_, z_, z_levels)
                    my_data['qs'][:, :, :, ctime] = interplevel(qs_, z_, z_levels)
                    my_data['qg'][:, :, :, ctime] = interplevel(qg_, z_, z_levels)

                except:
                    print('Warning no condensate data found')

                thetae_ = getvar(ncfile, "theta_e", units='K', timeidx=itime)
                theta_ = getvar(ncfile, "theta", units='K', timeidx=itime)
                rh_ = getvar(ncfile, "rh", timeidx=itime)
                try:
                    h_diabatic_ = getvar(ncfile, "H_DIABATIC", timeidx=itime)
                    qv_diabatic_ = getvar(ncfile, "QV_DIABATIC", timeidx=itime)
                    my_data['h_diabatic'][:, :, :, ctime] = interplevel(h_diabatic_, z_, z_levels)
                    my_data['qv_diabatic'][:, :, :, ctime] = interplevel(qv_diabatic_, z_, z_levels)
                except:
                    print('Warning no diabatic data found')
                my_data['p'][:, :, :, ctime] = interplevel(p_, z_, z_levels)
                my_data['u'][:, :, :, ctime] = interplevel(u_, z_, z_levels)
                my_data['v'][:, :, :, ctime] = interplevel(v_, z_, z_levels)
                my_data['w'][:, :, :, ctime] = interplevel(w_, z_, z_levels)
                my_data['t'][:, :, :, ctime] = interplevel(t_, z_, z_levels)
                my_data['qv'][:, :, :, ctime] = interplevel(qv_, z_, z_levels)
                my_data['thetae'][:, :, :, ctime] = interplevel(thetae_, z_, z_levels)
                my_data['theta'][:, :, :, ctime] = interplevel(theta_, z_, z_levels)
                my_data['rh'][:, :, :, ctime] = interplevel(rh_, z_, z_levels)


                my_data['times'].append(file_times[itime])

        if slice_type == 'h':
            for iz, my_z in enumerate(z_levels):
                my_data['z'][iz, :, :, :] = my_z

        if nt == 1:
            my_data['dt'] = 1.0
        else:
            my_data['dt'] = (my_data['times'][1] - my_data['times'][0]).total_seconds()
        my_data['dx'] = ncfile.DX
        my_data['dy'] = ncfile.DY
        my_data['file_name'] = file_name
        my_data['slice_type'] = slice_type
        my_data['slice_index'] = slice_index
        my_data['slice_width'] = slice_width
        my_data['slice_z'] = slice_z
        my_data['zres'] = zres
        my_data['ts'] = ts
        my_data['te'] = te

        my_data['p'] = np.float32(my_data['p'])
        my_data['u'] = np.float32(my_data['u'])
        my_data['v'] = np.float32(my_data['v'])
        my_data['w'] = np.float32(my_data['w'])
        my_data['ref'] = np.float32(my_data['ref'])
        my_data['t'] = np.float32(my_data['t'])
        my_data['qv'] = np.float32(my_data['qv'])
        my_data['qc'] = np.float32(my_data['qc'])
        my_data['qr'] = np.float32(my_data['qr'])
        my_data['qi'] = np.float32(my_data['qi'])
        my_data['qs'] = np.float32(my_data['qs'])
        my_data['qg'] = np.float32(my_data['qg'])
        my_data['z'] = np.float32(my_data['z'])
        my_data['thetae'] = np.float32(my_data['thetae'])
        my_data['theta'] = np.float32(my_data['theta'])
        my_data['rh'] = np.float32(my_data['rh'])
        my_data['h_diabatic'] = np.float32(my_data['h_diabatic'])
        my_data['qv_diabatic'] = np.float32(my_data['qv_diabatic'])

        if write_pkl:
            with open(file_name, 'wb') as handle:
                pickle.dump(my_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return my_data


def get_profile(data_path, xp=0, yp=0, tp=0):
    my_data = dict()

    my_data['file_list'] = glob.glob(data_path + '/wrfout_d01_*')  # Busco todos los wrfout en la carpeta indicada.
    my_data['file_list'].sort()

    ncfile = Dataset(my_data['file_list'][tp])
    my_data['p'] = to_np(getvar(ncfile, "pressure"))[:, yp, xp]
    my_data['z'] = to_np(getvar(ncfile, "height"))[:, yp, xp]
    my_data['t'] = to_np(getvar(ncfile, "temp", units='degC'))[:, yp, xp]
    my_data['td'] = to_np(getvar(ncfile, "td", units='degC'))[:, yp, xp]
    my_data['u'] = to_np(getvar(ncfile, "ua", units='ms-1'))[:, yp, xp]
    my_data['v'] = to_np(getvar(ncfile, "va", units='ms-1'))[:, yp, xp]
    my_data['xp'] = xp
    my_data['yp'] = yp
    my_data['tp'] = tp
    my_data['qv'] = to_np(getvar(ncfile, "QVAPOR"))[:, yp, xp]
    my_data['thetae'] = to_np(getvar(ncfile, "theta_e", units='K'))[:, yp, xp]
    my_data['theta'] = to_np(getvar(ncfile, "theta", units='K'))[:, yp, xp]

    return my_data


def get_profile_single_file(my_file, xp=0, yp=0, tp=0):
    my_data = dict()

    # my_data['file_list'] = glob.glob(data_path + '/wrfout_d01_*')   #Busco todos los wrfout en la carpeta indicada.
    # my_data['file_list'].sort()

    ncfile = Dataset(my_file)
    my_data['p'] = to_np(getvar(ncfile, "pressure", timeidx=tp))[:, yp, xp]
    my_data['z'] = to_np(getvar(ncfile, "height", timeidx=tp))[:, yp, xp]
    my_data['t'] = to_np(getvar(ncfile, "temp", units='degC', timeidx=tp))[:, yp, xp]
    my_data['td'] = to_np(getvar(ncfile, "td", units='degC', timeidx=tp))[:, yp, xp]
    my_data['u'] = to_np(getvar(ncfile, "ua", units='ms-1', timeidx=tp))[:, yp, xp]
    my_data['v'] = to_np(getvar(ncfile, "va", units='ms-1', timeidx=tp))[:, yp, xp]
    my_data['xp'] = xp
    my_data['yp'] = yp
    my_data['tp'] = tp
    my_data['qv'] = to_np(getvar(ncfile, "QVAPOR", timeidx=tp))[:, yp, xp]
    my_data['thetae'] = to_np(getvar(ncfile, "theta_e", units='K', timeidx=tp))[:, yp, xp]
    my_data['theta'] = to_np(getvar(ncfile, "theta", units='K', timeidx=tp))[:, yp, xp]

    return my_data


def get_ref_state(my_data, force=False):
    from scipy.interpolate import interp1d

    # Definimos el estado de referencia como el estado de las variables en el tiempo 0 para la esquina inferior izquierda del dominio que ingresa
    # El dominio que ingresa puede ser el dominio total o un corte.
    # No calculo el estado de referencia para los condensados (lo vamos a asumir siempre como 0).

    variables = ['p', 't', 'u', 'v', 'qv', 'theta', 'thetae']

    # Chequeamos si la cuenta ya se hizo.
    is_done = True
    for my_var in variables:
        if not my_var + '0' in my_data.keys():
            is_done = False

            # Si no se hizo la hacemos.
    if (not is_done) or force:
        print('Computing base state')

        # Para facilitar calculos posteriores tomamos el perfil de referencia y lo repetimos tantas veces como dimensiones tengan las variables.
        # Esto lo hacemos solo para el tiempo 0.

        for my_var in variables:

            if my_data['slice_type'] != 'h':
                var0 = my_data[my_var][:, 0, 0, 3]
                my_data[my_var + '0'] = np.zeros((my_data[my_var].shape))
                int_z = my_data['z'][:, 0, 0,
                        3]  # Uso el tiempo 1 en lugar del 0 ya que en el tiempo inicial por algun motivo el primer nivel de P esta mal.
                int_z[0] = int_z[0] - 1.0;
                int_z[-1] = int_z[-1] + 1.0
                interpolator = interp1d(int_z, var0, bounds_error=False, fill_value=np.nan)
                for ii in range(my_data[my_var].shape[1]):
                    for jj in range(my_data[my_var].shape[2]):
                        for kk in range(my_data[my_var].shape[3]):
                            my_data[my_var + '0'][:, ii, jj, kk] = interpolator(my_data['z'][:, ii, jj, kk])
            else:

                var0 = np.copy(my_data[my_var][:, 0, 0, 0])
                nrep = my_data[my_var].shape[1] * my_data[my_var].shape[2] * my_data[my_var].shape[3]
                tmp = np.tile(var0, nrep)
                tmp = np.reshape(tmp, (
                my_data[my_var].shape[1], my_data[my_var].shape[2], my_data[my_var].shape[3], var0.size))
                my_data[my_var + '0'] = tmp.transpose(3, 0, 1, 2)
    return my_data


def get_moment_equation(my_data, force=False, save=False):
    # Calculamos los terminos de la ecuacion de movimiento
    # Lista de variables que se agregan en esta funcion.
    variables = ['dwdt_loc', 'dwdt_adv', 'dwdt_bt', 'dwdt_bp', 'dwdt_bqv', 'dwdt_bqc', 'dwdt_pz', 'dudt_loc',
                 'dudt_adv', 'dudt_px', 'dvdt_loc', 'dvdt_adv', 'dvdt_py']

    # Chequeamos si la cuenta ya se hizo.
    is_done = True
    for my_var in variables:
        if not my_var in my_data.keys():
            is_done = False

    # Si no se hizo la hacemos.
    if (not is_done) or force:
        print('Computing momentum equation terms')

        my_data = get_ref_state(my_data, force=force)

        nz = my_data['p'].shape[0]
        nx = my_data['p'].shape[2]
        ny = my_data['p'].shape[1]
        nt = my_data['p'].shape[3]
        dt = my_data['dt']
        dx = my_data['dx']

        # Allocate memory
        for my_var in variables:
            my_data[my_var] = np.zeros((nz, ny, nx, nt)).astype('float32')

        # Calculo la derivada temporal local.
        if nt > 1:
            my_data['dwdt_loc'] = np.gradient(my_data['w'], axis=3) / my_data['dt']
            my_data['dudt_loc'] = np.gradient(my_data['u'], axis=3) / my_data['dt']
            my_data['dvdt_loc'] = np.gradient(my_data['v'], axis=3) / my_data['dt']

        for it in range(nt):
            u = my_data['u'][:, :, :, it]
            v = my_data['v'][:, :, :, it]
            w = my_data['w'][:, :, :, it]
            t = my_data['t'][:, :, :, it]
            t0 = my_data['t0'][:, :, :, it]
            tpert = t - t0
            p = my_data['p'][:, :, :, it]
            p0 = my_data['p0'][:, :, :, it]
            ppert = p - p0
            z = my_data['z'][:, :, :, it]
            qv = my_data['qv'][:, :, :, it]
            qv0 = my_data['qv0'][:, :, :, it]
            qvpert = qv - qv0
            qc = my_data['qc'][:, :, :, it] + my_data['qr'][:, :, :, it] + my_data['qi'][:, :, :, it] + my_data['qs'][:,
                                                                                                        :, :, it] + \
                 my_data['qg'][:, :, :, it]

            if my_data['slice_type'] == 'h':
                zlev = True
                dz = my_data['zres']
            else:
                zlev = False
                dz = None

            # Calculo el terminio advectivo.
            [my_data['dwdt_adv'][:, :, :, it], tmp, tmp, tmp] = advection(w, u, v, w, z, dx, dx, dz=dz, zlev=zlev)
            [my_data['dudt_adv'][:, :, :, it], tmp, tmp, tmp] = advection(u, u, v, w, z, dx, dx, dz=dz, zlev=zlev)
            [my_data['dvdt_adv'][:, :, :, it], tmp, tmp, tmp] = advection(v, u, v, w, z, dx, dx, dz=dz, zlev=zlev)
            # Calculo la densidad
            ro = p / (Rd * t * (1.0 + 0.608 * qv))
            roinv = 1.0 / ro
            # Calculo los terminos del empuje
            my_data['dwdt_bt'][:, :, :, it] = gg * tpert / t0
            my_data['dwdt_bp'][:, :, :, it] = - gg * ppert / p0
            my_data['dwdt_bqv'][:, :, :, it] = gg * 0.608 * qvpert
            my_data['dwdt_bqc'][:, :, :, it] = - gg * qc

            # Calculo la fuerza del gradiente de presion.

            [pz, py, px] = gradient_3d(ppert, z, dx, dx, dz=dz, zlev=zlev)
            my_data['dudt_px'][:, :, :, it] = -roinv * px
            my_data['dvdt_py'][:, :, :, it] = -roinv * py
            my_data['dwdt_pz'][:, :, :, it] = -roinv * pz

        if save:
            with open(my_data['file_name'], 'wb') as handle:
                pickle.dump(my_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return my_data


def get_ppert_equation(my_data, force=False, save=False):
    # Calculamos los terminos de la ecuacion de perturbaciones de presion.
    # Lista de variables que se agregan en esta funcion.
    variables = ['lp_spin_tot', 'lp_splat_tot', 'lp_spin_pert', 'lp_splat_pert', 'lp_b', 'lp_lineal', 'lp_ppert',
                 'bouy', 'ppert_nhy', 'ppert_hy', 'fppertx', 'fpperty', 'fppertz']
    #              'pp_spin_tot','pp_splat_tot','pp_spin_pert','pp_splat_pert','pp_b','pp_lineal','pp_tot']

    # Chequeamos si la cuenta ya se hizo.
    is_done = True
    for my_var in variables:
        if not my_var in my_data.keys():
            is_done = False

    # Si no se hizo la hacemos.
    if (not is_done) or force:
        print('Computing P\' equation terms')

        my_data = get_ref_state(my_data, force=force)

        nz = my_data['p'].shape[0]
        ny = my_data['p'].shape[1]
        nx = my_data['p'].shape[2]
        nt = my_data['p'].shape[3]
        dt = my_data['dt']
        dx = my_data['dx']

        # Allocate memory
        for my_var in variables:
            my_data[my_var] = np.zeros((nz, ny, nx, nt)).astype('float32')

        for it in range(nt):
            u = my_data['u'][:, :, :, it]
            v = my_data['v'][:, :, :, it]
            u0 = my_data['u0'][:, :, :, it]
            v0 = my_data['v0'][:, :, :, it]
            upert = u - u0
            vpert = v - v0
            w = my_data['w'][:, :, :, it]
            t = my_data['t'][:, :, :, it]
            t0 = my_data['t0'][:, :, :, it]
            tpert = t - t0
            p = my_data['p'][:, :, :, it]
            p0 = my_data['p0'][:, :, :, it]
            ppert = p - p0
            z = my_data['z'][:, :, :, it]
            qv = my_data['qv'][:, :, :, it]
            qv0 = my_data['qv0'][:, :, :, it]
            qvpert = qv - qv0
            qc = my_data['qc'][:, :, :, it] + my_data['qr'][:, :, :, it] + my_data['qi'][:, :, :, it] + my_data['qs'][:,
                                                                                                        :, :, it] + \
                 my_data['qg'][:, :, :, it]
            b = gg * tpert / t0 - gg * ppert / p0 + gg * 0.608 * qvpert - gg * qc
            ro0 = p0 / (Rd * t0 * (1.0 + 0.608 * qv0))
            ro = p / (Rd * t * (1.0 + 0.608 * qv))
            ropert = ro - ro0
            # Laplaciano de las perturbaciones de presion.
            my_data['lp_ppert'][:, :, :, it] = laplacian_3d(ppert, z, dx, dx)

            # Terminos dinamicos totales
            if my_data['slice_type'] == 'h':
                zlev = True
                dz = my_data['zres']
            else:
                zlev = False
                dz = None
            [uz, uy, ux] = gradient_3d(u, z, dx, dx, dz=dz, zlev=zlev)
            [vz, vy, vx] = gradient_3d(v, z, dx, dx, dz=dz, zlev=zlev)
            [wz, wy, wx] = gradient_3d(w, z, dx, dx, dz=dz, zlev=zlev)
            my_data['lp_spin_tot'][:, :, :, it] = 0.5 * (ro0) * ((wy - vz) ** 2 + (uz - wx) ** 2 + (vx - uy) ** 2)
            my_data['lp_splat_tot'][:, :, :, it] = -0.25 * (ro0) * (
                        (2.0 * ux) ** 2 + 2.0 * (uy + vx) ** 2 + 2.0 * (uz + wx) ** 2 + (2.0 * vy) ** 2 + 2.0 * (
                            vz + wy) ** 2 + (2.0 * wz) ** 2)
            # Terminos dinamicos perturbados
            [uz, uy, ux] = gradient_3d(upert, z, dx, dx, dz=dz, zlev=zlev)
            [vz, vy, vx] = gradient_3d(vpert, z, dx, dx, dz=dz, zlev=zlev)
            my_data['lp_spin_pert'][:, :, :, it] = 0.5 * (ro0) * ((wy - vz) ** 2 + (uz - wx) ** 2 + (vx - uy) ** 2)
            my_data['lp_splat_pert'][:, :, :, it] = -0.25 * (ro0) * (
                        (2.0 * ux) ** 2 + (uy + vx) ** 2 + (uz + wx) ** 2 + (vx + uy) ** 2 + (2.0 * vy) ** 2 + (
                            vz + wy) ** 2 + (wx + uz) ** 2 + (wy + vz) ** 2 + (2.0 * wz) ** 2)
            # Termino lineal.
            zz = np.gradient(z, axis=0)
            u0z = np.gradient(u0, axis=0) / zz
            v0z = np.gradient(v0, axis=0) / zz
            my_data['lp_lineal'][:, :, :, it] = -2.0 * (ro0) * (u0z * wx + v0z * wy)
            # Termino por empuje.
            [my_data['lp_b'][:, :, :, it], tmp, tmp] = gradient_3d(b * ro0, z, dx, dx, dz=dz, zlev=zlev)
            my_data['bouy'][:, :, :, it] = b

            # Calculo la perturbacion de presion hidrostatica.
            # Asumo que las perturbaciones en el tope del dominio son puramente hidrostaticas
            # En el tope del dominio predominan las ondas de gravedad donde la perturbacion tiene una componente hidrostatica importante.
            # Tomo como condicion inicial de la perturbacion hidrostatica la perturbacion de presion en el tope del dominio e integro de ahi para abajo.
            # Primero tomo el primer nivel que no sea nan en pprima (ese va a ser la perturbacion en el tope)
            ppert_top = np.zeros(ppert.shape[1:])
            for ii in range(ppert.shape[1]):
                for jj in range(ppert.shape[2]):
                    tmp_ppert = ppert[:, ii, jj]
                    tmp_ppert = tmp_ppert[~np.isnan(tmp_ppert)]
                    ppert_top[ii, jj] = tmp_ppert[-1]
            my_data['ppert_hy'][:, :, :, it] = np.tile(ppert_top, (ppert.shape[0], 1, 1))
            ropert_dz = 0.5 * np.diff(z, axis=0) * (ropert[0:-1, :, :] + ropert[1:, :, :]) * gg
            ropert_dz[np.isnan(ropert_dz)] = 0.0
            my_data['ppert_hy'][0:-1, :, :, it] = my_data['ppert_hy'][0:-1, :, :, it] + np.flip(
                np.cumsum(np.flip(ropert_dz), axis=0))
            my_data['ppert_nhy'][:, :, :, it] = ppert - my_data['ppert_hy'][:, :, :, it]

            # Calculo las componentes de la fuerza asociada a las perturbaciones de presion
            [pz, py, px] = gradient_3d(ppert, z, dx, dx, dz=dz, zlev=zlev)
            my_data['fppertx'][:, :, :, it] = -(1.0 / ro0) * px
            my_data['fpperty'][:, :, :, it] = -(1.0 / ro0) * py
            my_data['fppertz'][:, :, :, it] = -(1.0 / ro0) * pz

            # Solve the eliptic equation to get pressure perturbations associated to different forcings.
            # max_error=1.0e-4
            # if my_data['slice_type'] == 'vy' or my_data['slice_type'] == 'vx' :   #Veritcal slice (apply Dirichlet at the top)
            #    nws=[True,True,True]
            #    nwe=[False,True,True]
            #    ds=[False,False,False]
            #    de=[True,False,False]
            # if my_data['slice_type'] == 'h'  :                                    #Horizontal slice (apply Dirichlet at the lateral boundaries)
            #    nws=[True,False,False]
            #    nwe=[True,False,False]
            #    ds=[False,True,True]
            #    de=[False,True,True]

            # my_data['pp_spin_tot'][:,:,:,it] = solve_eliptic( my_data['lp_spin_tot'][:,:,:,it] , z , dx , dx , max_error = max_error , newmans=nws , newmane=nwe , dirichlets=ds , dirichlete=de )
            # my_data['pp_splat_tot'][:,:,:,it] = solve_eliptic( my_data['lp_splat_tot'][:,:,:,it] , z , dx , dx , max_error = max_error , newmans=nws , newmane=nwe , dirichlets=ds , dirichlete=de )
            # my_data['pp_spin_pert'][:,:,:,it] = solve_eliptic( my_data['lp_spin_pert'][:,:,:,it] , z , dx , dx , max_error = max_error , newmans=nws , newmane=nwe , dirichlets=ds , dirichlete=de )
            # my_data['pp_splat_pert'][:,:,:,it] = solve_eliptic( my_data['lp_splat_pert'][:,:,:,it] , z , dx , dx , max_error = max_error , newmans=nws , newmane=nwe , dirichlets=ds , dirichlete=de )
            # my_data['pp_b'][:,:,:,it] = solve_eliptic( my_data['lp_b'][:,:,:,it] , z , dx , dx , max_error = max_error , newmans=nws , newmane=nwe , dirichlets=ds , dirichlete=de )
            # my_data['pp_lineal'][:,:,:,it] = solve_eliptic( my_data['lp_lineal'][:,:,:,it] , z , dx , dx , max_error = max_error , newmans=nws , newmane=nwe , dirichlets=ds , dirichlete=de )
            # #my_data['pp_tot'][:,:,:,it] = solve_eliptic( my_data['lp_tot'][:,:,:,it] , z , dx , dx , max_error = max_error , newmans=nws , newmane=nwe , dirichlets=ds , dirichlete=de )

        if save:
            with open(my_data['file_name'], 'wb') as handle:
                pickle.dump(my_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return my_data


def get_vorticity_equation(my_data, force=False, save=False):
    # Calculamos los terminos de la ecuacion de perturbaciones de presion.
    # Lista de variables que se agregan en esta funcion.
    variables = ['vortx', 'vorty', 'vortz', 'dvortxdt_loc', 'dvortydt_loc', 'dvortzdt_loc', 'advz_vortx', 'advz_vorty',
                 'advz_vortz', 'advh_vortx', 'advh_vorty', 'advh_vortz', 'tilt_vortx', 'tilt_vorty', 'tilt_vortz',
                 'str_vortx', 'str_vorty', 'str_vortz', 'bouy_vortx', 'bouy_vorty']

    # Chequeamos si la cuenta ya se hizo.
    is_done = True
    for my_var in variables:
        if not my_var in my_data.keys():
            is_done = False

    # Si no se hizo la hacemos.
    if (not is_done) or force:
        print('Computing vorticity equation terms')

        my_data = get_ref_state(my_data, force=force)

        nz = my_data['p'].shape[0]
        ny = my_data['p'].shape[1]
        nx = my_data['p'].shape[2]
        nt = my_data['p'].shape[3]
        dt = my_data['dt']
        dx = my_data['dx']

        # Allocate memory
        for my_var in variables:
            my_data[my_var] = np.zeros((nz, ny, nx, nt)).astype('float32')

        for it in range(nt):
            u = my_data['u'][:, :, :, it]
            v = my_data['v'][:, :, :, it]
            u0 = my_data['u0'][:, :, :, it]
            v0 = my_data['v0'][:, :, :, it]
            w = my_data['w'][:, :, :, it]
            t = my_data['t'][:, :, :, it]
            t0 = my_data['t0'][:, :, :, it]
            tpert = t - t0
            p = my_data['p'][:, :, :, it]
            p0 = my_data['p0'][:, :, :, it]
            ppert = p - p0
            z = my_data['z'][:, :, :, it]
            qv = my_data['qv'][:, :, :, it]
            qv0 = my_data['qv0'][:, :, :, it]
            qvpert = qv - qv0
            qc = my_data['qc'][:, :, :, it] + my_data['qr'][:, :, :, it] + my_data['qi'][:, :, :, it] + my_data['qs'][:,
                                                                                                        :, :, it] + \
                 my_data['qg'][:, :, :, it]
            b = gg * tpert / t0 - gg * ppert / p0 + gg * 0.608 * qvpert - gg * qc

            if my_data['slice_type'] == 'h':
                zlev = True
                dz = my_data['zres']
            else:
                zlev = False
                dz = None

            # Calculo la vorticidad
            [uz, uy, ux] = gradient_3d(u, z, dx, dx, dz=dz, zlev=zlev)
            [vz, vy, vx] = gradient_3d(v, z, dx, dx, dz=dz, zlev=zlev)
            [wz, wy, wx] = gradient_3d(w, z, dx, dx, dz=dz, zlev=zlev)

            vortx = (wy - vz)
            vorty = (uz - wx)
            vortz = (vx - uy)

            my_data['vortx'][:, :, :, it] = vortx
            my_data['vorty'][:, :, :, it] = vorty
            my_data['vortz'][:, :, :, it] = vortz

            # Calculo los terminos de tilting
            my_data['tilt_vortx'][:, :, :, it] = uy * vorty + uz * vortz
            my_data['tilt_vorty'][:, :, :, it] = vx * vortx + vz * vortz
            my_data['tilt_vortz'][:, :, :, it] = wx * vortx + wy * vorty

            # Calculo los terminos de stretching
            my_data['str_vortx'][:, :, :, it] = -vortx * (vy + wz)
            my_data['str_vorty'][:, :, :, it] = -vorty * (ux + wz)
            my_data['str_vortz'][:, :, :, it] = -vortz * (ux + vy)

            # Calculo los terminos advectivos
            [tmp, advx_vx, advy_vx, advz_vx] = advection(vortx, u, v, w, z, dx, dx, dz=dz, zlev=zlev)
            [tmp, advx_vy, advy_vy, advz_vy] = advection(vorty, u, v, w, z, dx, dx, dz=dz, zlev=zlev)
            [tmp, advx_vz, advy_vz, advz_vz] = advection(vortz, u, v, w, z, dx, dx, dz=dz, zlev=zlev)

            my_data['advz_vortx'][:, :, :, it] = advz_vx
            my_data['advz_vorty'][:, :, :, it] = advz_vy
            my_data['advz_vortz'][:, :, :, it] = advz_vz
            my_data['advh_vortx'][:, :, :, it] = advx_vx + advy_vx
            my_data['advh_vorty'][:, :, :, it] = advx_vy + advy_vy
            my_data['advh_vortz'][:, :, :, it] = advx_vz + advy_vz

            # Caluclo los terminos de empuje.
            [bz, by, bx] = gradient_3d(b, z, dx, dx)
            my_data['bouy_vortx'][:, :, :, it] = by
            my_data['bouy_vorty'][:, :, :, it] = bx

        # Calculo la derivada temporal local.
        if nt > 1:
            my_data['dvortxdt_loc'] = np.gradient(my_data['vortx'], axis=3) / my_data['dt']
            my_data['dvortydt_loc'] = np.gradient(my_data['vorty'], axis=3) / my_data['dt']
            my_data['dvortzdt_loc'] = np.gradient(my_data['vortz'], axis=3) / my_data['dt']

        if save:
            with open(my_data['file_name'], 'wb') as handle:
                pickle.dump(my_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return my_data


def get_termo_equation(my_data, force=False, save=False):
    # Calculamos los terminos de la ecuacion de movimiento
    # Lista de variables que se agregan en esta funcion.
    variables = ['dthetadt_loc', 'dthetadt_adv']

    # Chequeamos si la cuenta ya se hizo.
    is_done = True
    for my_var in variables:
        if not my_var in my_data.keys():
            is_done = False

    # Si no se hizo la hacemos.
    if (not is_done) or force:
        print('Computing termodinamic equation terms')

        my_data = get_ref_state(my_data, force=force)

        nz = my_data['p'].shape[0]
        ny = my_data['p'].shape[1]
        nx = my_data['p'].shape[2]
        nt = my_data['p'].shape[3]
        dt = my_data['dt']
        dx = my_data['dx']

        my_data['dthetadt_loc'] = np.zeros((nz, ny, nx, nt)).astype(
            'float32')  # Derivada local de theta con respecto al tiempo
        my_data['dthetadt_adv'] = np.zeros((nz, ny, nx, nt)).astype('float32')  # Termino advectivo.

        # Calculo la derivada temporal local.
        if nt > 1:
            my_data['dthetadt_loc'] = np.gradient(my_data['theta'], axis=3) / dt

        for it in range(nt):
            u = my_data['u'][:, :, :, it]
            v = my_data['v'][:, :, :, it]
            w = my_data['w'][:, :, :, it]
            theta = my_data['theta'][:, :, :, it]
            z = my_data['z'][:, :, :, it]

            if my_data['slice_type'] == 'h':
                zlev = True
                dz = my_data['zres']
            else:
                zlev = False
                dz = None

            # Calculo el terminio advectivo.
            [my_data['dthetadt_adv'][:, :, :, it], tmp, tmp, tmp] = advection(theta, u, v, w, z, dx, dx, dz=dz,
                                                                              zlev=zlev)

        if save:
            with open(my_data['file_name'], 'wb') as handle:
                pickle.dump(my_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return my_data


def get_water_equation(my_data, force=False, save=False):
    # Calculamos los terminos de la ecuacion de movimiento
    # Lista de variables que se agregan en esta funcion.
    variables = ['dqvdt_loc', 'dqvdt_adv', 'dqvdt_diabatic', 'dqdt_loc', 'dqcdt_adv', 'dqrdt_loc', 'dqrdt_adv',
                 'dqidt_loc', 'dqidt_adv', 'dqsdt_loc', 'dqsdt_adv', 'dqgdt_loc', 'dqgdt_adv']

    # Chequeamos si la cuenta ya se hizo.
    is_done = True
    for my_var in variables:
        if not my_var in my_data.keys():
            is_done = False

    # Si no se hizo la hacemos.
    if (not is_done) or force:
        print('Computing water equation terms')

        my_data = get_ref_state(my_data, force=force)

        nz = my_data['p'].shape[0]
        ny = my_data['p'].shape[1]
        nx = my_data['p'].shape[2]
        nt = my_data['p'].shape[3]
        dt = my_data['dt']
        dx = my_data['dx']

        for my_var in variables:
            my_data[my_var] = np.zeros((nz, ny, nx, nt)).astype('float32')

        # Calculo la derivada temporal local.
        if nt > 1:
            my_data['dqvdt_loc'] = np.gradient(my_data['qv'], axis=3) / dt
            my_data['dqcdt_loc'] = np.gradient(my_data['qc'], axis=3) / dt
            my_data['dqrdt_loc'] = np.gradient(my_data['qr'], axis=3) / dt
            my_data['dqidt_loc'] = np.gradient(my_data['qi'], axis=3) / dt
            my_data['dqsdt_loc'] = np.gradient(my_data['qs'], axis=3) / dt
            my_data['dqgdt_loc'] = np.gradient(my_data['qg'], axis=3) / dt

        for it in range(nt):
            u = my_data['u'][:, :, :, it]
            v = my_data['v'][:, :, :, it]
            w = my_data['w'][:, :, :, it]
            z = my_data['z'][:, :, :, it]
            qv = my_data['qv'][:, :, :, it]
            qc = my_data['qc'][:, :, :, it]
            qr = my_data['qr'][:, :, :, it]
            qi = my_data['qi'][:, :, :, it]
            qs = my_data['qs'][:, :, :, it]
            qg = my_data['qg'][:, :, :, it]

            if my_data['slice_type'] == 'h':
                zlev = True
                dz = my_data['zres']
            else:
                zlev = False
                dz = None

            # Calculo el terminio advectivo.
            [my_data['dqvdt_adv'][:, :, :, it], tmp, tmp, tmp] = advection(qv, u, v, w, z, dx, dx, dz=dz, zlev=zlev)
            [my_data['dqcdt_adv'][:, :, :, it], tmp, tmp, tmp] = advection(qc, u, v, w, z, dx, dx, dz=dz, zlev=zlev)
            [my_data['dqrdt_adv'][:, :, :, it], tmp, tmp, tmp] = advection(qr, u, v, w, z, dx, dx, dz=dz, zlev=zlev)
            [my_data['dqidt_adv'][:, :, :, it], tmp, tmp, tmp] = advection(qi, u, v, w, z, dx, dx, dz=dz, zlev=zlev)
            [my_data['dqsdt_adv'][:, :, :, it], tmp, tmp, tmp] = advection(qs, u, v, w, z, dx, dx, dz=dz, zlev=zlev)
            [my_data['dqgdt_adv'][:, :, :, it], tmp, tmp, tmp] = advection(qg, u, v, w, z, dx, dx, dz=dz, zlev=zlev)

        if save:
            with open(my_data['file_name'], 'wb') as handle:
                pickle.dump(my_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return my_data


def gradient_3d(var, z, dy, dx, dz=None, zlev=False):
    # El gradiente es a z constante (aun cuando las niveles verticales no lo sean, eg sigma o p)
    # Calculo primero la derivada vertical que es mas compleja porque los puntos no estan equiespaciados.

    if not zlev:  # Irregular grid spacing in z
        varz = np.zeros(var.shape).astype('float32')
        for ii in range(var.shape[1]):
            for jj in range(var.shape[2]):
                varz[:, ii, jj] = np.gradient(var[:, ii, jj], z[:, ii, jj], axis=0)

        vary = np.gradient(var, dy, axis=1)
        varx = np.gradient(var, dx, axis=2)

        zy = np.gradient(z, dy, axis=1)
        zx = np.gradient(z, dx, axis=2)

        varx = (varx - zx * varz)
        vary = (vary - zy * varz)

    else:  # Constant z levels
        [varz, vary, varx] = np.gradient(var)
        varz = varz / dz
        varx = varx / dx
        vary = vary / dy

    return varz, vary, varx


def laplacian_3d(var, z, dy, dx):
    # Calculo el laplaciano como la divergencia del gradiente y uso la funcion
    # gradient de numpy. Las derivadas horizontales estan corregidas para tener
    # en cuanta los cambios de z con x y con y en superficies como p y sigma.

    [varz, vary, varx] = gradient_3d(var, z, dy, dx)
    [varzz, varzy, varzx] = gradient_3d(varz, z, dy, dx)
    zx = np.gradient(z, dx, axis=2)
    zy = np.gradient(z, dy, axis=1)
    # Uso que varxz = varzx y que varyz = varzy
    varxx = np.gradient(varx, dx, axis=2) - zx * varzx
    varyy = np.gradient(vary, dy, axis=1) - zy * varzy

    return varzz + varyy + varxx


def advection(var, u, v, w, z, dy, dx, dz=None, zlev=False):
    [varz, vary, varx] = gradient_3d(var, z, dy, dx, dz=dz, zlev=zlev)
    advx = -u * varx
    advy = -v * vary
    advz = -w * varz

    return advx + advy + advz, advx, advy, advz


# def solve_eliptic( F , z , dy , dx , max_error = 1.e-3 , newmans=[True,True,True] , newmane=[True,True,True] , dirichlets=[False,False,False] , dirichlete=[False,False,False] ) :
# #Vamos a resolver la ecuacion eliptica con las siguientes condiciones:
# #F es el forzante dentro del dominio.
# #La condicion de borde superior es fija (Dirichlet) y es p=0
# #La condicion de borde a los lados y en la base es de Newman y es
# #dp_estrelladx=0.
# #Usamos el metodo de relajacion (Jacobi) para resolver la ecuacion.

#     sol_error=1.0
#     sol_ = np.zeros( np.shape( F ) )
#     [nz,ny,nx] = np.shape(F)
#     max_error=5.0e-4

#     dx_ = 1.0*dx
#     dy_ = 1.0*dy
#     dz_ = ( z[2:,1:-1,1:-1] - z[1:-1,1:-1,1:-1] ) * 1.0
#     factor = ( 2.0/ (dx_**2) + 2.0/ (dy_**2) + 2.0/(dz_**2) )


#     while (sol_error > max_error ) :
#        #Primero resuelvo en el interior del dominio.
#        lap_ = laplacian_3d( sol_ , z , dy , dx )
#        #Actualizo el valor de p
#        sol_[1:-1,1:-1,1:-1]=sol_[1:-1,1:-1,1:-1] - 1.6*(1.0/factor) * ( F[1:-1,1:-1,1:-1] - ( lap_[1:-1,1:-1,1:-1] ) )
#        sol_error=(np.abs(F-lap_)[1:-1,1:-1,1:-1]).max()
#        print(sol_error)

#        if newmans[0] :
#           sol_[0,:,:] = sol_[1,:,:]
#        if newmane[0] :
#           sol_[nz-1,:,:] = sol_[nz-2,:,:]
#        if newmans[1] :
#           sol_[:,0,:] = sol_[:,1,:]
#        if newmane[1] :
#           sol_[:,ny-1,:] = sol_[:,ny-2,:]
#        if newmans[2] :
#           sol_[:,:,0] = sol_[:,:,1]
#        if newmans[2] :
#           sol_[:,:,nx-1] = sol_[:,:,nx-2]

#        if dirichlets[0] :
#           sol_[0,:,:] = 0.0
#        if dirichlete[0] :
#           sol_[nz-1,:,:] = 0.0
#        if dirichlets[1] :
#           sol_[:,0,:] = 0.0
#        if dirichlete[1] :
#           sol_[:,ny-1,:] = 0.0
#        if dirichlets[2] :
#           sol_[:,:,0] = 0.0
#        if dirichlete[2] :
#           sol_[:,:,nx-1] = 0.0


#     return sol_

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
         max_u = conf['shear_depth_v'] * conf['shear_strength_v']
         input_sounding['v'] = input_sounding['height'] * conf['shear_strength_v']
         input_sounding['v'][ input_sounding['height'] > conf['shear_depth_v'] ] = max_v

      if conf['shear_type'] == 'Curved' :
         #Maximum turning in wind direction (is a function of how much shear is explained by the curved part of the hodograph)
         #The maximum turning is pi (half a circle), if half the shear is associated with the curved part then a quarter circle hodograph is obtained.
         curved_shear_theta = conf['curved_shear_per'] * np.pi
         #Curved hodograp up to certain level, then linear shear.

         #Curved hodograp part
         if conf['curved_shear_per'] > 0.0 :
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

      #print( input_sounding['u'] , input_sounding['v'] ) 


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
       #factor = 1.0 + 0.01 * conf['stability_factor'] *  input_sounding['height'] / ( 10000.0 ) 
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


def sensitivity_function_single_file(my_file, sf_type='near_surface_wind_percentile', sf_percentile=99, inilev=-1,
                                     endlev=-1, initime=-1, endtime=-1):
    # This function computes sensitivity metrics taking the full 4D model output (for different variables) and computing a scalar function of it that
    # represent some aspect of the sensitivity of the output to the model parameters.

    # Leo el primer archivo de la lista de donde voy a tomar el perfil de referencia.
    ncfile = Dataset(my_file)
    tmp = to_np(getvar(ncfile, "P"))
    nx = tmp.shape[2]  # Order is Z=0 , Y=1 , X=0
    ny = tmp.shape[1]
    nz = tmp.shape[0]
    raw_times = extract_times(ncfile, None, meta=False)
    # Convert numpy.datetime64 into datetime format.
    file_times = [dt.datetime.strptime(np.datetime_as_string(tmpdate, timezone='UTC')[0:19], '%Y-%m-%dT%H:%M:%S') for
                  tmpdate in raw_times]
    nt = len(file_times)
    if initime == -1:
        initime = 0
    if endtime == -1:
        endtime = nt - 1
    if inilev == -1:
        inilev = 0
    if endlev == -1:
        endlev = nz - 1

    if (sf_type == 'near_surface_wind_percentile'):
        # This scalar function consists on a percentile of the wind speed distribution (over time and space) for a given layer (for example the first 10 model levels).
        # Using a high percentile such as 99 % this function will focus on the most extreme values. However the use of a percentile is more stable and less sensitive
        # to noise than the global maximum.
        # WARNING: This parameter is affected by the mean wind speed which is assumed to be 0 in all simulations (to keep the storm within the domain as long as possible).
        # this can lead to unrealistic results since the actual wind speed at the surface will be given by the wind speed relative to the storm + the storm velocity (which
        # in this experients is artificially kept near 0). A possible correction would be to assume that at the initial time the wind speed at the surface is 0 and correct
        # the maximum wind speed by this factor (taking the sign into account).

        ncfile = Dataset(my_file)
        wspeed = np.zeros((nz, ny, nx, nt))
        for itime in range(nt):
            wspeed[:, :, :, itime] = np.sqrt(to_np(getvar(ncfile, "ua", units='ms-1', timeidx=itime)) ** 2 + to_np(
                getvar(ncfile, "va", units='ms-1', timeidx=itime)) ** 2)

        # Compute the requested percentile of the near surface wind.
        sf = np.percentile(np.ndarray.flatten(wspeed[inilev:endlev, :, :, initime:endtime]), sf_percentile)

    elif (sf_type == 'upward_vertical_velocity_percentile'):
        # This scalar function consists on a percentile of the upward vertical velocity. As a measure of how intense are the resulting convective cells.
        # Using a high percentile such as 99 % this function will focus on the most extreme values. However the use of a percentile is more stable and less sensitive
        # to noise than the global maximum.

        ncfile = Dataset(my_file)
        wa = np.zeros((nz, ny, nx, nt))
        for itime in range(nt):
            wa[:, :, :, itime] = to_np(getvar(ncfile, "wa", units='ms-1', timeidx=itime))

        # Compute the requested percentile of the near surface wind.
        wa = wa[inilev:endlev, :, :, initime:endtime]
        sf = np.percentile(wa[wa > 0], sf_percentile)


    elif (sf_type == 'mean_precipitation'):
        # This scalar function consists on a percentile of the upward vertical velocity. As a measure of how intense are the resulting convective cells.
        # Using a high percentile such as 99 % this function will focus on the most extreme values. However the use of a percentile is more stable and less sensitive
        # to noise than the global maximum.

        ncfile = Dataset(my_file)
        sf = np.mean(to_np(getvar(ncfile, "RAINNC", timeidx=endtime)) - to_np(getvar(ncfile, "RAINNC",
                                                                                     timeidx=initime)))  # Domain-wise and time averaged precipitation (mm , equivalent to liters / m^2 )

    else:
        print('Error: Not recognized sensitivity function')

    return sf


def edit_namelist(file_path, conf):
    import os

    # Create temp file
    # ftmp = open( file_path + '.tmp' )
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

    # Copy the file permissions from the old file to the new file
    os.system('rm -f ' + file_path)
    os.system('mv ' + file_path + '.tmp ' + file_path)


def run_wrf(conf):
    #############################################################################################
    # PREPARATORY STEPS
    #############################################################################################
    print('start wrf nunning from wrf module')
    # Create output directory
    expoutpath = conf['datapath'] + '/' + conf['expname'] + '/' + '{:03d}'.format( conf['run_num'] )
    os.makedirs(expoutpath, exist_ok=True)
    # Copy configuration file to the output directory
    with open(expoutpath + '/conf_' + '{:03d}'.format(conf['run_num']) + '.pkl', 'wb') as handle:
        pickle.dump(conf, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Create temporary directory
    # We need a temporary directory to run wrf model
    # we create a random one assuming that several experiments will be running at the same time.

    tmp_dir = conf['datapath'] + '/tmp/tmpdir_' + conf['expname'] + '{:03d}'.format( conf['run_num'] )
    os.makedirs(tmp_dir, exist_ok=True)

    #############################################################################################
    # STEP 1 - GENERATE THE STORM "ENVIRONMENT"
    #############################################################################################

    if conf['run_model']:
        # Copy wrf executables and data to the temp directory.
        os.system('cp ' + conf['modelpath'] + '/* ' + tmp_dir + '/')

        # Create the terrain
        if conf.get('terr_generate'):

            from terraingen import generate_terrain, terrain_to_string

            terrParams = {}
            for key in conf.keys():
                if key.startswith('terr_'):
                    terrParams[key[len('terr_'):]] = conf[key]
            terrain = generate_terrain(**terrParams)
            terrain_file = open(tmp_dir + '/terrain.txt', "w")
            terrain_file.write(terrain_to_string(terrain))
            terrain_file.close()

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

    #############################################################################################
    # STEP 4 - OBTAIN SIMULATION METRICS
    #############################################################################################
    '''
    if conf['run_model'] :

        sensitivity_score = sensitivity_function_single_file(outdatafile, sf_type=conf['sf_type'],
                                                             sf_percentile=conf['sf_percentile'],
                                                             inilev=conf['sf_inilev'], endlev=conf['sf_endlev']
                                                             , initime=conf['sf_initime'], endtime=conf['sf_endtime'])

        print('The sensitivity score is : ', sensitivity_score)
    '''
    #############################################################################################
    # STEP 5 (OPTIONAL) - PLOT THE OUTPUT
    #############################################################################################

    if conf['plot_exp']:
        # Import modules for ploting
        import wrf_plot as wrfp
        # Create a directory for the figures.
        plot_path = expoutpath + '/figs_' + '{:03d}'.format(conf['run_num']) + '/'
        os.makedirs(plot_path, exist_ok=True)

        # Plot time evolution of "reflectivity"
        wrf_data_2 = get_data_vslice_single_file( outdatafile , slice_width=1 , slice_z=1000.0 , zres=1 , slice_index=1 , slice_type='h', t_start=0 ,
                                                 t_end=0 , write_pkl=False )
        if 'ref_vertvel' in conf['plot_types']:
            wrfp.plot_ref_vertvel(wrf_data_2, plot_path, scale_factor=1.0, arrow_scale_factor=5.0)

        if 'horwind_theta' in conf['plot_types']:
            wrfp.plot_horwind_theta(wrf_data_2, plot_path, scale_factor=1.0, arrow_scale_factor=5.0)

        if 'init_sounding' in conf['plot_types']:
            # Plot storm "environment"
            wrf_data = get_profile_single_file(outdatafile, xp=0, yp=0, tp=0)  # Read the data from netcdf file.
            wrfp.plot_sounding(wrf_data, plot_path, show=False)  # Plot the storm "environment"

    #return sensitivity_score

