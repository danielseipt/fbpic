# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)
It defines a set of utility functions for particle initialization.
"""
import numpy as np
import warnings
from scipy.constants import c

def boost_and_select( x, y, z, ux, uy, uz, w, comm, boost ):
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

    boost : float or None
        The Lorentz factor of the boosted frame
    """
    inv_gamma = 1./np.sqrt( 1. + ux**2 + uy**2 + uz**2 )
    # Convert the particles to the boosted-frame
    if boost is not None:
        x, y, z, ux, uy, uz, inv_gamma = boost.boost_particle_arrays(
                                        x, y, z, ux, uy, uz, inv_gamma )

    # Select the particles that are in the local subdomain
    zmin, zmax = comm.get_zmin_zmax(
        local=True, with_damp=False, with_guard=False, rank=comm.rank )
    selected = (z >= zmin) & (z < zmax)
    x = x[selected]
    y = y[selected]
    z = z[selected]
    ux = ux[selected]
    uy = uy[selected]
    uz = uz[selected]
    w = w[selected]
    inv_gamma = inv_gamma[selected]

    return( x, y, z, ux, uy, uz, inv_gamma, w )


def generate_evenly_spaced(dens_func, v_m, v_th, n_per_cell, comm, zmin, zmax):
    """
    Generate evenly-spaced particles, according to the density function
    `dens_func`, and with the momenta given by the `v_m/v_th` arguments.

    n_per_cell: dictionary with keys 'r', 'z', 'theta'
    """
    # Use bounds of local domain, if zmin and zmax are not provided
    if zmin is None:
        zmin, _ = comm.get_zmin_zmax( local=True, with_damp=False,
                                    with_guard=False, rank=comm.rank )
    if zmax is None:
        zmax, _ = comm.get_zmin_zmax( local=True, with_damp=False,
                                    with_guard=False, rank=comm.rank )
    Npz = n_per_cell['z'] * int( round( (zmax-zmin)/comm.dz ) )
    Npr = n_per_cell['r'] * int( round( comm.rmax/comm.dr ) )
    Nptheta = n_per_cell['theta']

    # Generate the particles and eliminate the ones that have zero weight ;
    # infer the number of particles Ntot
    if Npz*Npr*Nptheta > 0:
        # Get the 1d arrays of evenly-spaced positions for the particles
        dz = (zmax-zmin)*1./Npz
        z_reg =  zmin + dz*( np.arange(Npz) + 0.5 )
        dr = comm.rmax*1./Npr
        r_reg =  rmin + dr*( np.arange(Npr) + 0.5 )
        dtheta = 2*np.pi/Nptheta
        theta_reg = dtheta * np.arange(Nptheta)

        # Get the corresponding particles positions
        # (copy=True is important here, since it allows to
        # change the angles individually)
        zp, rp, thetap = np.meshgrid( z_reg, r_reg, theta_reg,
                                    copy=True, indexing='ij' )
        # Prevent the particles from being aligned along any direction
        unalign_angles( thetap, Npz, Npr, method='random' )
        # Flatten them (This performs a memory copy)
        r = rp.flatten()
        x = r * np.cos( thetap.flatten() )
        y = r * np.sin( thetap.flatten() )
        z = zp.flatten()
        # Get the weights (i.e. charge of each macroparticle), which
        # are equal to the density times the volume r d\theta dr dz
        w = r * dtheta*dr*dz * dens_func( z, r )

        # Select the particles that have a non-zero weight
        selected = (w > 0)
        if np.any(w < 0):
            warnings.warn(
            'The specified particle density returned negative densities.\n'
            'No particles were generated in areas of negative density.\n'
            'Please check the validity of the `dens_func`.')

        # Infer the number of particles and select them
        Ntot = int(selected.sum())
        x = x[ selected ]
        y = y[ selected ]
        z = z[ selected ]
        w = w[ selected ]
        # Initialize the corresponding momenta
        ux = np.random.normal( v_m[0]/c, v_th[0]/c, Ntot)
        uy = np.random.normal( v_m[1]/c, v_th[1]/c, Ntot)
        uz = np.random.normal( v_m[2]/c, v_th[2]/c, Ntot)
        inv_gamma = 1./np.sqrt( 1 + ux**2 + uy**2 + uz**2 )
        # Return the particle arrays
        return( x, y, z, ux, uy, uz, inv_gamma, w )
    else:
        # No particles are initialized ; the arrays are still created
        return( np.empty(0), np.empty(0), np.empty(0), np.empty(0),
                np.empty(0), np.empty(0), np.empty(0), np.empty(0) )


def unalign_angles( thetap, Npz, Npr, method='irrational' ) :
    """
    Shift the angles so that the particles are
    not all aligned along the arms of a star transversely

    The fact that the particles are all aligned can produce
    numerical artefacts, especially if the polarization of the laser
    is aligned with this direction.

    Here, for each position in r and z, we add the *same*
    shift for all the Nptheta particles that are at this position.
    (This preserves the fact that certain modes are 0 initially.)
    How this shift varies from one position to another depends on
    the method chosen.

    Parameters
    ----------
    thetap : 3darray of floats
        An array of shape (Npr, Npz, Nptheta) containing the angular
        positions of the particles, and which is modified by this function.

    Npz, Npr : ints
        The number of macroparticles along the z and r directions

    method : string
        Either 'random' or 'irrational'
    """
    # Determine the angle shift
    if method == 'random' :
        angle_shift = 2*np.pi*np.random.rand(Npz, Npr)
    elif method == 'irrational' :
        # Subrandom sequence, by adding irrational number (sqrt(2) and sqrt(3))
        # This ensures that the sequence does not wrap around and induce
        # correlations
        shiftr = np.sqrt(2)*np.arange(Npr)
        shiftz = np.sqrt(3)*np.arange(Npz)
        angle_shift = 2*np.pi*( shiftz[:,np.newaxis] + shiftr[np.newaxis,:] )
        angle_shift = np.mod( angle_shift, 2*np.pi )
    else :
        raise ValueError(
      "method must be either 'random' or 'irrational' but is %s" %method )

    # Add the angle shift to thetap
    # np.newaxis ensures that the angles that are at the same positions
    # in r and z have the same shift
    thetap[:,:,:] = thetap[:,:,:] + angle_shift[:,:, np.newaxis]
