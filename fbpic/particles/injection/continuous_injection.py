# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)
It defines a class for continuous particle injection with a moving window.
"""
import numpy as np
from scipy.constants import c
from .particle_layouts import GriddedLayout

class ContinuousInjector( object ):
    """
    Class that stores a number of attributes that are needed for
    continuous injection by a moving window.
    """

    def __init__(self, distribution, layout, comm, boost ):
        """
        Initialize continuous injection

        Parameters
        ----------
        #TODO add parameters
        """
        # Register the spacing between particles in z
        if isinstance( layout, GriddedLayout ):
            self.dz_particles = comm.dz/layout.n_macroparticle_per_cell['z']
        else:
            self.dz_particles = 0
        # Register the boost
        # (necessary when generating particles directly in the boosted frame)
        self.boost = boost

        # Register variables that define the positions
        # where the plasma is injected.
        v_m = np.array( getattr(distribution, 'directed_velocity', [0,0,0]) )
        self.v_end_plasma = v_m[2] / np.sqrt(1 + np.dot(v_m,v_m)/c**2)
        # These variables are set by `initialize_injection_positions`
        self.nz_inject = None
        self.z_inject = None
        self.z_end_plasma = None


    def initialize_injection_positions( self, comm, v_moving_window,
                                        species_z, dt ):
        """
        Initialize the positions that keep track of the injection of particles.
        This is automatically called at the beginning of `step`.

        Parameters
        ----------
        comm: a BoundaryCommunicator object
            Contains information about grid MPI decomposition
        v_moving_window: float (in m/s)
            The speed of the moving window
        species_z: 1darray of float (in m)
            (One element per macroparticle)
            Used in order to infer the position of the end of the plasma
        dt: float (in s)
            Timestep of the simulation
        """
        # The injection position is only ini
        if comm.rank != comm.size-1:
            return
        # Initialize the injection position only if it has not be initialized
        if self.z_inject is not None:
            return

        # Initialize plasma *ahead* of the right *physical*
        # boundary of the box so that, after `exchange_period` iterations
        # (without adding new plasma), there will still be plasma
        # inside the physical domain. ( +3 takes into account that
        # 3 more cells need to be filled w.r.t the left edge of the
        # physical box such that the last cell inside the box is
        # always correct for 1st and 3rd order shape factor
        # particles after the moving window shifted by exchange_period cells.)
        _, zmax_global_domain = comm.get_zmin_zmax( local=False,
                                    with_damp=False, with_guard=False )
        self.z_inject = zmax_global_domain + 3*comm.dz + \
                comm.exchange_period*dt*(v_moving_window-self.v_end_plasma)
        self.nz_inject = 0
        # Try to detect the position of the end of the plasma:
        # Find the maximal position of the continously-injected particles
        if len( species_z ) > 0:
            # Add half of the spacing between particles (the
            # injection function itself will add a half-spacing again)
            self.z_end_plasma = species_z.max() + 0.5*self.dz_particles
        else:
            # Default value for empty species
            self.z_end_plasma = zmax_global_domain

        # Check that the particle spacing has been properly calculated
        if self.dz_particles is None:
            raise ValueError(
                'The simulation uses continuous injection of particles, \n'
                'but was unable to calculate the spacing between particles.\n'
                'This may be because you used the `Particles` API directly.\n'
                'In this case, please pass the argument `dz_particles` \n'
                'initializing the `Particles` object.')

    def reset_injection_positions( self ):
        """
        Reset the variables that keep track of continuous injection to `None`
        This is typically called when restarting a simulation from a checkpoint
        """
        self.nz_inject = None
        self.z_inject = None
        self.z_end_plasma = None

    def increment_injection_positions( self, v_moving_window, duration ):
        """
        Update the positions between which the new particles will be generated,
        the next time when `generate_particles` is called.
        This function is automatically called when the moving window moves.

        Parameters
        ----------
        v_moving_window: float (in m/s)
            The speed of the moving window

        duration: float (in seconds)
            The duration since the last time that the moving window moved.
        """
        # Move the injection position
        self.z_inject += v_moving_window * duration
        # Take into account the motion of the end of the plasma
        self.z_end_plasma += self.v_end_plasma * duration

        # Increment the number of particle to add along z
        nz_new = int( (self.z_inject - self.z_end_plasma)/self.dz_particles )
        self.nz_inject += nz_new
        # Increment the virtual position of the end of the plasma
        # (When `generate_particles` is called, then the plasma
        # is injected between z_end_plasma - nz_inject*dz_particles
        # and z_end_plasma, and afterwards nz_inject is set to 0.)
        self.z_end_plasma += nz_new * self.dz_particles


    def generate_particles( self, distribution, layout, comm, time ):
        """
        Generate new particles at the right end of the plasma
        (i.e. between z_end_plasma - nz_inject*dz and z_end_plasma)

        TODO: update parameters

        Parameters
        ----------
        time: float (in second)
            The current physical time of the simulation
        """
        # Create new particle cells
        # Determine the positions between which new particles will be created
        zmax = self.z_end_plasma
        zmin = self.z_end_plasma - self.nz_inject*self.dz_particles
        # Create the particles
        x, y, z, ux, uy, uz, inv_gamma, w = distribution.generate_particles(
            layout, comm, self.boost, time, zmin, zmax)

        # Reset the number of particle cells to be created
        self.nz_inject = 0

        return( x, y, z, ux, uy, uz, inv_gamma, w )
