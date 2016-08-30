import qdune
import os
import numpy as np


__all__ = ['dotest']

def dotest():
    # File name formatting string things
    datadir = qdune.__path__[0]+os.sep+'testdata'

    # input arguments
    save = False
    flow_azimuth = 158
    dx = 0.1
    dy = 0.25
    istart = 1800
    iend = 2300
    Lmax = (iend-istart)*dx/3
    ref_xyz = 0

    qdat = qdune.q(datadir, save, flow_azimuth, dx, dy, istart, iend, Lmax)

    flux2d = qdat.Vc * qdat.Hc * 0.5
    volflux = np.nansum(flux2d,1) * qdat.dy
    massflux = volflux * 2650 * 0.65

    print('Average Flux = {0} kg/s'.format(np.nanmean(massflux)))
