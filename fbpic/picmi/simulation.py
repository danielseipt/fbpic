# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)

It defines the picmi Simulation and CylindricalGrid interface
"""
from fbpic.main import Simulation as FBPICSimulation
from fbpic.smoothing import BinomialSmoother # Exposed directly, unmodified

class CylindricalGrid(object):

    def __init__(self, nr, rmin, rmax, nz, zmin, zmax, n_azimuthal_modes,
                bc_zmin='periodic', bc_zmax='periodic', bc_rmax='reflective',
                moving_window_velocity=0. ):
        """
        TODO
        """
        self.nr = nr
        self.rmin = rmin
        self.rmax = rmax
        self.nz = nz
        self.zmin = zmin
        self.zmax = zmax
        if bc_rmax != 'reflective':
            raise ValueError(
            "In fbpic, `bc_rmax` can only be 'reflective' - not %s." %bc_rmax)
        self.bc_rmax = bc_rmax
        if bc_zmin not in ['periodic', 'open']:
            raise ValueError("In fbpic, `bc_zmin` can only be "
                "'periodic' or 'open' - not %s." %bc_zmin)
        self.bc_zmin = bc_zmin
        if bc_zmax != bc_zmax:
            raise ValueError(
            "In fbpic, `bc_zmax` should be identical to `bc_zmin`.")
        self.bc_zmax = bc_zmax
        self.n_azimuthal_modes = n_azimuthal_modes
        self.moving_window_velocity = moving_window_velocity


class ElectromagneticSolver(object):

    def __init__(self, grid, cfl=None, dt=None, source_smoother=None ):
        """
        TODO
        """
        self.grid = grid
        self.cfl = cfl
        self.dt = dt
        self.source_smoother = source_smoother
        # TODO Check that either cfl or dt is passed ; automatically calculate dt


class Simulation(FBPICSimulation):

    def __init__(self, solver, dt=None, **kw ):
        """
        TODO
        """
        grid = solver.grid
        if dt is None:
            dt = solver.dt
        else:
            # TODO Print warning message
            pass

        # Initialize the object using the standard FBPIC API
        FBPICSimulation.__init__(self,
            Nz=grid.nz, zmin=grid.zmin, zmax=grid.zmax,
            Nr=grid.nr, rmin=grid.rmin, rmax=grid.rmax,
            dt=dt, boundaries=grid.bc_zmin,
            smoother=solver.source_smoother, **kw )
        # TODO: Not currently implemented in PICMI:
        # galilean, gamma_boost, verbose_level, use_all_mpi_ranks, use_cuda

        # Remove the default empty electron species in fbpic
        self.ptcl = []
