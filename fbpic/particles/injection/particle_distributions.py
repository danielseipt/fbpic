# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)
It defines a set of physical particle distributions.
"""
import warnings
import numexpr
import numpy as np
from scipy.constants import c
from .particle_layouts import PseudoRandomLayout, GriddedLayout
from .utilities import boost_and_select, generate_evenly_spaced

class DistributionFromArrays( object ):

    def __init__( self, x, y, z, vx, vy, vz, w ):
        """
        TODO

        vx, vy, vz: velocities (gamma*v)

        Numpy Arrays of positions in the lab frame
        """
        self.x = x
        self.y = y
        self.z = z
        self.ux = vx/c
        self.uy = vy/c
        self.uz = vz/c
        self.w = w

    def generate_particles( self, layout, comm, boost, **kw ):
        """Returns x, y, z, ux, uy, uz, inv_gamma, w as arrays"""

        if layout is not None:
            warnings.warn('The parameter `layout` is ignored '
                          'when using `DistributionFromArrays`.')
        return boost_and_select( self.x, self.y, self.z,
                self.ux, self.uy, self.uz, self.w, comm, boost )


class GaussianBunchDistribution(object):

    def __init__(self, n_physical_particles, rms_bunch_size,
        rms_velocity=[0.,0.,0.], centroid_position=[0.,0.,0.],
        centroid_velocity=[0.,0.,0.], velocity_divergence=[0.,0.,0.]):
        """
        TODO
        velocities are expressed in terms of gamma*v
        """
        self.n_physical_particles = n_physical_particles
        self.rms_bunch_size = rms_bunch_size
        self.rms_velocity = rms_velocity
        self.centroid_position = centroid_position
        self.centroid_velocity = centroid_velocity
        self.velocity_divergence = velocity_divergence

    def generate_particles( self, layout, comm, boost, **kw ):
        """Returns x, y, z, ux, uy, uz, inv_gamma, w as arrays"""

        # Randomly generate the bunch
        if isinstance(layout, PseudoRandomLayout):
            N = layout.n_macroparticles
            # Get Gaussian particle distribution in x,y,z
            position = self.centroid_position
            bunch_size = self.rms_bunch_size
            x = np.random.normal( position[0], bunch_size[0], N )
            y = np.random.normal( position[1], bunch_size[1], N )
            z = np.random.normal( position[2], bunch_size[2], N )
            # Get Gaussian particle distribution in ux, uy, uz
            v = self.centroid_velocity
            dv = self.rms_velocity
            ux = np.random.norma( v[0]/c, dv[0]/c, N )
            uy = np.random.norma( v[1]/c, dv[1]/c, N )
            uz = np.random.norma( v[2]/c, dv[2]/c, N )
            # TODO: implement velocity divergence
            w = self.n_physical_particles/N * np.ones(N)

            # Return the particle arrays
            return boost_and_select(x, y, z, ux, uy, uz, w, comm, boost)

        else:
            raise ValueError('`layout` of type %s is not implemented '
                            'for GaussianBunchDistribution.' %type(layout) )


class PythonFunctionDistribution( object ):

    def __init__(self, density_function,
                 lower_bound = [None,None,None],
                 upper_bound = [None,None,None],
                 rms_velocity_spread = [0.,0.,0.],
                 directed_velocity = [0.,0.,0.],
                 fill_in = False, **kw):
        self.density_function = density_function
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.rms_velocity_spread = rms_velocity_spread
        self.directed_velocity = directed_velocity
        self.fill_in = fill_in

    def generate_particles( self, layout, comm, boost, zmin=None, zmax=None ):
        """Returns x, y, z, ux, uy, uz, inv_gamma, w as arrays"""
        if isinstance(layout, GriddedLayout):

            return generate_evenly_spaced( self.density_function,
                    self.directed_velocity, self.rms_velocity_spread,
                    layout.n_macroparticle_per_cell, comm, zmin, zmax )

            # TODO: Handle boost of input length, etc. ; same as WarpX
            # (However, this would perturb the current workings of FBPIC,
            # which expects the user to boost the lengths.)
            # TODO: Scale the density of the particles
            # TODO: Boost the velocity of the particles

        else:
            raise ValueError('`layout` of type %s is not implemented '
                            'for AnalyticDistribution.' %type(layout) )


class AnalyticDistribution( PythonFunctionDistribution ):

    def __init__(self, density_expression,
                 lower_bound = [None,None,None],
                 upper_bound = [None,None,None],
                 rms_velocity_spread = [0.,0.,0.],
                 directed_velocity = [0.,0.,0.],
                 fill_in = False, **kw):

        # Parse the string using numexpr
        def density_function( x, y, z ):
            return numexpr.evaluate(density_expression)
        # Initialize a standard PythonFunctionDistribution
        PythonFunctionDistribution.__init__(self,
            density_function=density_function,
            lower_bound=lower_bound, upper_bound=upper_bound,
            rms_velocity_spread=rms_velocity_spread,
            directed_velocity=directed_velocity, fill_in=fill_in )
