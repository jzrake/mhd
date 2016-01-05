import os
import numpy as np
import h5py


def load_tseries(rundir):
    data = np.loadtxt(os.path.join(rundir, 'mhd.dat'))
    return dict(
        t  = data[:,0],
        Eu = data[:,1],
        Hu = data[:,2],
        Mu = data[:,3],
        Lu = data[:,4],
        Eb = data[:,5],
        Hb = data[:,6],
        Mb = data[:,7],
        Lb = data[:,8])


def load_user_params(rundir):
    chkpt = os.path.join(rundir, 'chkpt.0000.h5')
    userf = os.path.join(rundir, 'user.h5')
    if os.path.isfile(chkpt):
        h5f = h5py.File(chkpt, 'r')
    else:
        h5f = h5py.File(userf, 'r')
    user = dict([(n, h5f['user'][n][0]) for n in h5f['user'].dtype.names])
    h5f.close()
    return user
