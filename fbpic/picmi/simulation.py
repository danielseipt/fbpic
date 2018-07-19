# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)

It defines the picmi Simulation and CylindricalGrid interface
"""
# Import relevant fbpic object
from fbpic.main import Simulation as FBPICSimulation
from fbpic.lpa_utils.laser import add_laser_pulse
from fbpic.lpa_utils.laser import GaussianLaser as FBPIC_GaussianLaser
# Import redefined picmi classes
from .laser import LaserAntenna, LaserDirect, GaussianLaser

class Simulation(FBPICSimulation):

    def __init__(self, solver, dt=None, gamma_boost=None, **kw ):
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
            gamma_boost=gamma_boost,
            smoother=solver.source_smoother, **kw )
        # TODO: Not currently implemented in PICMI:
        # galilean, gamma_boost, verbose_level, use_all_mpi_ranks, use_cuda

        # Remove the default empty electron species in fbpic
        self.ptcl = []


    def add_laser(self, laser, injection_method=None ):
        """
        TODO

        If injection_method is None, the default is laser antenna, in fbpic
        """
        # Convert injection method to corresponding fbpic arguments
        if isinstance(injection_method, LaserAntenna):
            method = 'antenna'
            z0_antenna = injection_method.position[2]
            v_antenna = injection_method.velocity[2]
        elif (injection_method is None) or \
                isinstance(injection_method, LaserDirect):
            method = 'direct'
            z0_antenna = None
            v_antenna = None
        else:
            raise ValueError("`injection_method` of type %s is not implemented."
                                %type(injection_method))

        # Convert laser to corresponding fbpic profile
        if isinstance(laser, GaussianLaser):
            profile = FBPIC_GaussianLaser(
                a0=laser.a0, waist=laser.waist, tau=laser.duration,
                z0=laser.centroid_position[2], zf=laser.focal_position[2],
                theta_pol=laser.polarization_angle, lambda0=laser.wavelength,
                propagation_direction=laser.propagation_direction[2] )
        else:
            raise ValueError(
                "`laser` of type %s is not recognized." %type(laser))

        # Call the proper fbpic function
        # (use the default `gamma_boost` of the simulation)
        add_laser_pulse(self, profile, method=method, z0_antenna=z0_antenna,
                        v_antenna=v_antenna, gamma_boost=self.gamma_boost)
