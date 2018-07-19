# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)
It defines a set of physical particle distributions.
"""
import warnings
import numpy as np
from scipy.constants import c
from fbpic.lpa_utils.boost import BoostConverter
from .particle_layouts import PseudoRandomLayout

# Particle distributions
# ----------------------

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

    def generate_particles( self, layout, comm, gamma_boost ):
        """Returns x, y, z, ux, uy, uz, inv_gamma, w as arrays"""

        if layout is not None:
            warnings.warn('The parameter `layout` is ignored '
                          'when using `DistributionFromArrays`.')
        return boost_and_select( self.x, self.y, self.z,
                self.ux, self.uy, self.uz, self.w, comm, gamma_boost )


class GaussianBunchDistribution(object):

    def __init__(n_physical_particles, rms_bunch_size,
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

    def generate_particles( self, layout, comm, gamma_boost ):
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
            return boost_and_select( x, y, z, ux, uy, uz, comm, gamma_boost )

        else:
            raise ValueError('`layout` of type %s is not implemented '
                            'for GaussianBunchDistribution.' %type(layout) )




# Helper functions
# ----------------

def boost_and_select( x, y, z, ux, uy, uz, w, comm, gamma_boost ):
    """
    TODO

    Parameters
    ----------
    x, y, z: 1d arrays of length (N_macroparticles,)
        The positions of the particles in x, y, z in meters

    ux, uy, uz: 1d arrays of length (N_macroparticles,)
        The dimensionless momenta of the particles in each direction

    w: 1d array of length (N_macroparticles,)
        The weight of the particles, i.e. the number of physical particles
        that each macroparticle corresponds to.

    gamma_boost : float or None
        The Lorentz factor of the boosted frame
    """
    inv_gamma = 1./np.sqrt( 1. + ux**2 + uy**2 + uz**2 )
    # Convert the particles to the boosted-frame
    if gamma_boost is not None:
        boost = BoostConverter(gamma_boost)
        x, y, z, ux, uy, uz, inv_gamma = boost.boost_particle_arrays(
                                        x, y, z, ux, uy, uz, inv_gamma )

    # Select the particles that are in the local subdomain
    zmin, zmax = comm.get_zmin_zmax(
        local=True, with_damp=False, with_guard=False, rank=sim.comm.rank )
    selected = (z >= zmin) & (z < zmax)
    x = x[selected]
    y = y[selected]
    z = z[selected]
    ux = ux[selected]
    uy = uy[selected]
    uz = uz[selected]
    w = w[selected]
    inv_gamma = inv_gamma[selected]

    return( x, y, z, ux, uy, uyz, inv_gamma, w )
