"""
Suite of tools for gridding multibeam bathymetric surveys and calculating
bed form flux.


version: 1.0
"""

import os, glob
import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import datetime as dt

class Raw(object):

    def __init__(self, datadir):
        """ Instantiates a multibeam raw data object """
        self.datadir = datadir

        # Load file info and verify data validity
        info = pd.read_csv(os.path.join(self.datadir, 'info.csv'))
        nt = len(info)

        self.t = [dt.datetime.strptime(str(info.ix[i,'DATE'])
                                      + str(info.ix[i,'TIME']),
                                      '%Y%m%d%H%M')
                  for i in range(nt)]

        self.xyzdir = os.path.join(datadir, 'xyz')
        self.xyznames = info.ix[:,'FILENAME']

        self.xyzpaths = [os.path.join(self.xyzdir, self.xyznames[i])
            for i in range(nt)]

        xyzglobpath = os.path.join(self.datadir, 'xyz', '*')
        assert (set(self.xyzpaths).issubset(set(glob.glob(xyzglobpath)))), \
        "FILENAME listed in 'info.csv' not found in datadir"

    # ================================================================
    def _rotate_xy(self, xyz, angle):
        """ Rotates every point in an xyz array about the origin by 'angle' """
        theta = np.deg2rad(angle)
        x = xyz[:,0].copy()
        y = xyz[:,1].copy()
        xyz[:,0] = x * np.cos(theta) - y * np.sin(theta)
        xyz[:,1] = x * np.sin(theta) + y * np.cos(theta)
        xyz[:,2] = xyz[:,2]
        return xyz

    # ================================================================
    def _generate_grid(self, xyz, dx, dy):
        """Gengenerate_grid(xyz, dx, dy)erate location data based on input XYZ data"""

        origin = np.amin(xyz,0)
        extent = np.amax(xyz,0)-origin
        ncells = (np.amax(xyz,0)-origin)//[dx,dy,1]

        # Account for remainder
        origin += [(extent[0] % dx) / 2, (extent[1] % dy) / 2, 0]

        xbnds = np.linspace(0, ncells[0] * dx, ncells[0] + 1)
        ybnds = np.linspace(0, ncells[1] * dy, ncells[1] + 1)

        return origin, xbnds, ybnds, extent[2]

    # ================================================================
    def _gridloc(self, nddata, bnds, axis):
        """Return grid indices of every point in xyz"""

        assert(axis < nddata.ndim), "axis > ndim"

        nddata = nddata[nddata[:,axis].argsort()]
        loc = np.searchsorted(bnds, nddata[:,axis]) - 1
        return nddata, loc

    # ================================================================
    def _smooth_profile(self, x, y, xbnds, usepts, nanrad):
        """ Smooth a single profile of 1D unstructured data """

        xz = np.hstack((x[:,np.newaxis],y[:,np.newaxis]))
        xz, xloc = self._gridloc(xz, xbnds, 0)

        dx = np.diff(xbnds)
        xnew = xbnds[:-1] + 0.5 * dx

        nanmask = np.array([np.sum(np.abs(xloc - ix+1) <= nanrad) == 0 \
                            for ix in range(len(xnew))])

        frac = usepts/xz.shape[0]

        with np.errstate(invalid='ignore'):
            try:
                w = sm.nonparametric.lowess(xz[:,1], xz[:,0], frac=frac,
                                            delta=np.nanmean(dx), it=0)
                znew = np.interp(xnew, w[:,0], w[:,1])
            except:
                znew = np.nan

        znew[nanmask] = np.nan
        return xnew, znew

    # ================================================================
    def _grid_xyz(self, xyz, xbnds, ybnds, zmax, usepts, nanrad):
        nx = len(xbnds) - 1
        ny = len(ybnds) - 1

        xyz, yloc = self._gridloc(xyz, ybnds, 1)

        Z = np.empty((ny, nx), 'float')
        for i in range(ny):

            xz = xyz[np.ix_(yloc==i, [0,2])]
            if len(xz) < usepts:
                Z[i,:] = np.nan
            else:
                with np.errstate(invalid='ignore'):
                    _, Z[i] = self._smooth_profile(xz[:,0], xz[:,1], xbnds,
                                             usepts=usepts,
                                             nanrad=nanrad)

        Z[np.ma.fix_invalid(Z, fill_value = np.inf) > zmax] = np.nan
        Z[np.ma.fix_invalid(Z, fill_value = 0) < 0] = np.nan

        return Z

    # ================================================================
    def plot_timestep(self, survey, rot_deg=0, markersize=5):
        """Load a single xyz file and plot datapoints colored by elevation.

        Intended to assist in determining flow azimuth.
        """
        xyz = self._rotate_xy(np.loadtxt(self.xyzpaths[survey]), rot_deg)

        plt.clf()
        fig1 = plt.figure(1)
        ax = fig1.add_subplot(111)

        sc = ax.scatter(xyz[:,0], xyz[:,1], c = xyz[:,2],
                        marker = 'o',
                        s = markersize, lw = 0,
                        cmap = 'viridis')

        fig1.colorbar(sc)
        ax.axis('equal')
        ax.set_title('Raw data for survey {0}'.format(survey))
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        plt.show()
        return ax


    # ================================================================
    def plot_smooth(self, survey, xsect, rot_deg, dx=1, dy=1, usepts=20, nanrad=2):
        """Plot raw data and smoothed transect for QC.

        Use this function to test gridding parameters.
        """

        xyz = self._rotate_xy(np.loadtxt(self.xyzpaths[survey]), rot_deg)
        origin, xbnds, ybnds, zmax = self._generate_grid(xyz, dx, dy)
        xyz -= origin
        xyz, yloc = self._gridloc(xyz, ybnds, 1)
        xz = xyz[np.ix_(yloc==xsect, [0,2])]

        xnew, znew = self._smooth_profile(xz[:,0], xz[:,1], xbnds, usepts, nanrad)

        plt.clf()
        fig1 = plt.figure(1)
        ax = fig1.add_subplot(111)

        ax.scatter(xz[:,0], xz[:,1])
        ax.plot(xnew, znew, 'r')

        ax.set_title('Transect Smoothing')
        ax.set_xlabel('Streamwise (x) coordinate')
        ax.set_ylabel('Z coordinate')
        plt.show()
        return ax

    def plot_gridded(self, survey, rot_deg, dx=1, dy=1, usepts=20, nanrad=2):
        """ Plot a single gridded survey for QC """

        # load data
        xyz = self._rotate_xy(np.loadtxt(self.xyzpaths[survey]), rot_deg)
        origin, xbnds, ybnds, zmax = self._generate_grid(xyz, dx, dy)
        xyz -= origin
        xyz, yloc = self._gridloc(xyz, ybnds, 1)

        # grid data
        zarr = self._grid_xyz(xyz, xbnds, ybnds, zmax, usepts=usepts, nanrad=2)
        xarr, yarr = np.meshgrid(xbnds[:-1], ybnds[:-1])

        # plot
        fig1 = plt.figure(1)
        ax = fig1.add_subplot(111)
        ax.pcolormesh(xarr, yarr, zarr, cmap = 'viridis', vmin = 0, vmax = zmax)
        ax.colorbar()
        ax.axis('tight')
        plt.show()

        return ax

# ============================================================================

class Raster(Raw):

    def __init__(self, datadir, save, ref_xyz, flow_azimuth, dx, dy,
                 usepts, nanrad):
        """ Return a gridded data object """
        Raw.__init__(self, datadir)
        self.ref_xyz = ref_xyz
        self.flow_azimuth = flow_azimuth
        self.dx = dx
        self.dy = dy

        # Generate Grid based on reference survey
        refdat = self._rotate_xy(np.loadtxt(self.xyzpaths[ref_xyz]), self.flow_azimuth)
        self.origin, self.xbnds, self.ybnds, self.zmax = \
                self._generate_grid(refdat, self.dx, self.dy)

        nt = len(self.t)
        nx = len(self.xbnds) - 1
        ny = len(self.ybnds) - 1

        # Make raster directory if it doesn't exist
        self.rasterdir = os.path.join(datadir, 'raster')

        if save == True and not os.path.exists(self.rasterdir):
            os.mkdir(self.rasterdir)

        # Generate list of raster pathnames
        rnames = ['r{0}.npy'.format(os.path.splitext(self.xyznames[i])[0]) \
                    for i in range(nt)]
        self.rpaths = [os.path.join(self.rasterdir, rnames[i]) for i in range(nt)]
        self.completed = glob.glob(os.path.join(self.rasterdir, '*'))

        # For every survey
        self.Zarr = np.empty((nt, ny, nx), 'float')
        for i in range(nt):

            # if the raster file already exists:
            if self.rpaths[i] in self.completed:
                # Load completed raster
                self.Zarr[i] = np.load(self.rpaths[i])

                print('\rVerified {0}/{1}'.format(i+1, nt), end="\r")
            else:
                # Load point clouds
                xyz = self._rotate_xy(np.loadtxt(self.xyzpaths[i]), flow_azimuth)
                xyz -= self.origin
                xyz, yloc = self._gridloc(xyz, self.ybnds, 1)

                # grid pointcloud
                zarr = self._grid_xyz(xyz, self.xbnds, self.ybnds,
                                              self.zmax, usepts=usepts,
                                              nanrad=nanrad)

                # save raster and store in Zarr
                if save == True:
                    np.save(self.rpaths[i], zarr)
                elif save == False:
                    pass
                else:
                    raise Exception('save must be boolean')

                self.Zarr[i] = zarr

                # record progress
                self.completed.append(self.rpaths[i])
                print('\rCompleted {0}/{1}'.format(i+1, nt), end = "\r")


    # =========================================================
    def _nan_helper(self, y):
        """ Helper to handle indices and logical indices of NaNs."""

        return np.isnan(y), lambda z: z.nonzero()[0]

    # =========================================================
    def _remove_outliers(self, data, nsigma):
        """ Iteratively removes outliers from data defined by nsigma """
        while (np.abs(np.ma.fix_invalid(data)-np.nanmean(data)) > nsigma *
        np.nanstd(data)).any():

            data[np.where(np.abs(np.ma.fix_invalid(data)-np.nanmean(data)) >
            nsigma * np.nanstd(data))] = np.nan

        return data

    # =========================================================
    def _div0(self, a, b , val = 0):
        """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
        with np.errstate(divide='ignore', invalid='ignore'):
            c = np.true_divide( a, b )
            c[ ~ np.isfinite( c )] = val  # -inf inf NaN
        return c

    # =========================================================
    def _xcorrf(self, profile1, profile2, dx):
        """ Find displacement and correlation coefficient of two bed form profiles"""
        corrf = np.correlate(profile2, profile1, mode = 'same') \
                /np.sum(profile1**2)

        if np.isnan(corrf).any():
            displ = np.nan
            corr = 0
        else:
            displ = (np.where(corrf == np.max(corrf))[0][0] - len(corrf)//2)*dx
            corr = np.max(corrf)

        return displ, corr

    # =========================================================
    def _lowpass(self, signal, dx, Lmax):
        """Return filtered signal according to maximum wavelength specification.

        Assumes signal is already zero mean and nans interpolated
        """
        W = np.fft.rfftfreq(len(signal), dx)
        f_signal = np.fft.rfft(signal)
        filtered_signal = f_signal.copy()
        filtered_signal[np.where(self._div0(1,W) > Lmax)] = 0+0j
        return np.fft.irfft(filtered_signal, n=len(signal)), f_signal

    # =========================================================
    def _detrend_xsect(self, xsect, dx, Lmax, mindat):
        """Return filtered transects based on maximum wavelength specification

        This function handles fourier filtering for bed form profiles with missing
        data.
        """
        nans,x = self._nan_helper(xsect)

        if sum(nans) > len(xsect)*(1-mindat):
            xsect[:] = np.nan
            f_signal = np.zeros(len(xsect)//2 + 1, 'complex')

        else:
            xsect[nans] = np.interp(x(nans), x(~nans), xsect[~nans])
            xsect -= np.nanmean(xsect)

            # filter signal
            xsect[:], f_signal = self._lowpass(xsect, dx, Lmax)

        xsect[nans] = np.nan
        return xsect, f_signal

    # =========================================================
    def _detrend_all(self, Z, dx, Lmax, datmin_grid):
        """Filter all transects based on maximum wavelength specification"""
        [nt, ny, nx] = Z.shape
        z = np.empty((nt, ny, nx), 'float')
        z_f = np.empty((nt, ny, nx//2 + 1), 'complex')

        for t in range(nt):
            for y in range(ny):
                z[t,y], z_f[t,y] = self._detrend_xsect(Z[t,y].copy(), dx, Lmax, datmin_grid)

        return z, z_f

    # =========================================================
    def _calc_Lc(self, signal, dx):
        """ Calculate characteristic lengthscale of bed form profile """

        Wwin = np.fft.rfftfreq(len(signal), dx)
        f_signal = np.fft.rfft(signal)
        amplitude = np.abs(f_signal) * dx
        power = amplitude ** 2
        fc = np.sum(power[1:]*Wwin[1:])/np.sum(power[1:])

        return 1/fc

    # =========================================================
    def _calc_Hc(self, signal):
        """ Calculate characteristic height scale of bed form profile"""

        return 2.8 * np.nanstd(signal)

    # =========================================================
    def _bf_geom(self, Zwin, dx, datmin_geom = 0.25, nsigma = 2.5, maxslope = 0.15):
        """ Calculate nt by ny bedform geometry array

        Replace unrealistic geometries with nan according to
        the following protocol:

        -remove Hc and Lc values that are greater than nsigma standard
        deviations from the mean value of timestep

        -remove values where Hc/Lc > maxslope

        -mask values that are only valid for Hc or Lc, so that Hc and Lc
        are valid simultaneously everywhere

        """
        [nt, ny, _] = Zwin.shape

        Hc = np.empty((nt, ny), 'float')
        Lc = np.empty((nt, ny), 'float')

        for t in range(nt):
            for y in range(ny):
                signal = Zwin[t,y]
                nans, x = self._nan_helper(signal)
                if sum(nans)/len(nans) < 1-datmin_geom:
                    Hc[t,y] = self._calc_Hc(signal)
                    signal[nans] = np.interp(x(nans), x(~nans), signal[~nans])
                    signal = signal - np.nanmean(signal)
                    Lc[t,y] = self._calc_Lc(signal, dx)
                else:
                    Hc[t,y] = np.nan
                    Lc[t,y] = np.nan
        # Remove outlier values
        for i in range(nt):
            Hc[i] = self._remove_outliers(Hc[i], nsigma)
            Lc[i] = self._remove_outliers(Lc[i], nsigma)

        # Remove unrealistic geometries
        ixbad = np.where(np.ma.fix_invalid(Hc/Lc) > maxslope)
        Hc[ixbad] = np.nan
        Lc[ixbad] = np.nan

        # Mask so Hc and Lc are always both valid
        nanmask = np.ma.mask_or(np.isnan(Hc), np.isnan(Lc))
        Hc[nanmask] = np.nan
        Lc[nanmask] = np.nan

        return Hc, Lc


    # =========================================================
    def _dx_mat(self, z, dx, mincorr):
        """ Return displacement matrix. """
        [nt, ny, _] = z.shape
        disp = np.empty((nt, nt, ny), 'float')
        for t1 in range(nt):
            for t2 in range(nt):
                for y in range(ny):
                    d, corr = self._xcorrf(z[t1,y], z[t2, y], dx)

                    if corr > mincorr:
                        disp[t1,t2,y] = d
                    else:
                        disp[t1, t2, y] = np.nan

            print("\rCorrelated t1 = {0}/{1}".format(t1+1, nt), end = "\r")

        return disp


    # =========================================================
    def _dt_mat(self, t):
        """ Return duration matrix. """
        nt = len(t)
        delta_t = np.empty((nt, nt), 'float')
        for t1 in range(nt):
            for t2 in range(nt):
                delta_t[t1,t2] = (t[t2] - t[t1]).total_seconds()

        return delta_t


    # =========================================================
    def _clean_dx_mat(self, disp, deltat, Vmin, Vmax, Lc, dfracmin, dfracmax):
        """ replace unrealistic velocities with nan according to the following
        conditions:

        - Vmax > Vc > Vmin
        - dfracmax > disp/Lc > dfracmin
        """
        [nt, _, ny] = disp.shape
        for t1 in range(nt):
            for t2 in range(nt):
                for y in range(ny):
                    d_x = disp[t1, t2, y]
                    d_t = deltat[t1,t2]

                    lc = Lc[t1,y]

                    ## comparison values
                    disp_max = np.min([lc * dfracmax,
                                      Vmax * abs(d_t)])
                    disp_min = np.max([lc * dfracmin,
                                     Vmin * abs(d_t)])

                    # test if valid
                    valid = disp_min < abs(d_x) < disp_max

                    # fix if not valid
                    if not valid:
                        disp[t1, t2, y] = np.nan
            print("\rCleaned {0}/{1}".format(t1+1, nt), end = "\r")

        return disp

    def _vregress(self, disp, deltat, r2min):
        """ perfoms velocity regression on valid displacements """
        [nt, _, ny] = disp.shape
        Vc = np.empty((nt, ny), 'float')
        r2 = np.empty_like(Vc)
        for t1 in range(nt):
            for y in range(ny):
                d_x = disp[t1,:,y].flatten()
                d_t = deltat[t1,:].flatten()

                ival = ~np.isnan(d_x)

                d_x = d_x[ival]
                d_t = d_t[ival]

                d_t = d_t[:,np.newaxis]

                if len(d_x) < 4:
                    Vc[t1, y] = np.nan
                else:
                    Vc[t1,y], resid, _, _ = np.linalg.lstsq(d_t, d_x)
                    r2[t1,y] = 1 - resid / (len(d_x) * np.var(d_x))
            print('\rVelocity Calculated {0}/{1}'.format(t1+1, nt), end = '\r')
        Vc[np.where(r2<r2min)] = np.nan

        return Vc

    def _bf_vel(self, z, t, dx, mincorr, Vmin, Vmax, Lc, dfracmin, dfracmax,
                minR2):
        """ Return nt by ny velocity array """
        displacement = self._dx_mat(z, dx, mincorr)
        duration = self._dt_mat(t)
        displacement = self._clean_dx_mat(displacement, duration,
                                     Vmin, Vmax, Lc, dfracmin, dfracmax)

        Vc = self._vregress(displacement, duration, minR2)
        return Vc

# ===================================================================

class FluxStats(Raster):

    def __init__(self, datadir, save, flow_azimuth, dx, dy, istart, iend, Lmax,
                 ref_xyz, datmin_grid, usepts, nanrad, datmin_geom, nsigma,
                 maxslope, mincorr, Vmin, Vmax, dfracmin, dfracmax, minR2):

        # Gridding parameters
        Raster.__init__(self, datadir, save, ref_xyz, flow_azimuth, dx, dy,
                             usepts, nanrad)
        self.istart = istart
        self.iend = iend
        self.Lmax = Lmax
        self.datmin_grid = datmin_grid

        # Geometric calculation parameters
        self.datmin_geom = datmin_geom
        self.nsigma = nsigma
        self.maxslope = maxslope

        # Velocity parameters
        self.mincorr = mincorr
        self.Vmin = Vmin
        self.Vmax = Vmax
        self.dfracmin = dfracmin
        self.dfracmax = dfracmax
        self.minR2 = minR2


        self.Zwin = self.Zarr[:,:,istart:iend]
        self.z, self.z_f = self._detrend_all(self.Zwin, dx, Lmax, datmin_grid)

        [self.Hc, self.Lc] = self._bf_geom(self.z, self.dx, self.datmin_geom,
                                           self.nsigma, self.maxslope)

        self.Vc = self._bf_vel(self.z, self.t, self.dx, self.mincorr,
                               self.Vmin, self.Vmax, self.Lc, self.dfracmin,
                               self.dfracmax, self.minR2)


def q(datadir, save, flow_azimuth, dx, dy, istart, iend, Lmax, ref_xyz=0,
      datmin_grid=0.25, usepts=20, nanrad=2, datmin_geom=0.8, nsigma=2.5,
      maxslope=0.15, mincorr=0.85, Vmin=0.1/3600, Vmax=3/3600, dfracmin=0.001,
      dfracmax=0.2, minR2=0.8):

    return FluxStats(datadir, save, flow_azimuth, dx, dy, istart, iend, Lmax,
                    ref_xyz, datmin_grid, usepts, nanrad, datmin_geom, nsigma,
                    maxslope, mincorr, Vmin, Vmax, dfracmin, dfracmax, minR2)
