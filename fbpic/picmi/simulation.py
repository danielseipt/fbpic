# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)

It defines the picmi Simulation and CylindricalGrid interface
"""
from fbpic.main import Simulation as FBPICSimulation
from fbpic.lpa_utils.laser import add_laser_pulse

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


    def add_laser(self, laser, injection ):

        if isinstance(injection, )

        add_laser_pulse(self, laser, )
