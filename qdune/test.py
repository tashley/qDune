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

    qobj = qdune.q(datadir, save, flow_azimuth, dx, dy, istart, iend, Lmax)

    flux2d = qobj.Vc * qobj.Hc * 0.5 * 0.65
    volflux = np.nansum(flux2d,1) * qobj.dy
    massflux = volflux * 2650

    correct_result = np.array([0.73350211, 0.83383889, 0.81234798, 0.51997741])
    assert not np.isnan(massflux).any(), "NaNs found in flux output"
    assert np.allclose(massflux, correct_result), "Flux calculation incorrect"
    print('Test Successful.', end='\t\t\t\t\n')
    return qobj
