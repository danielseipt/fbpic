# Copyright 2017, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)
It defines the structure that damps the fields in the guard cells.
"""
import numpy as np
try:
    from fbpic.cuda_utils import cuda_tpb_bpg_2d, cuda
    cuda_installed = cuda.is_available()
except ImportError:
    cuda_installed = False

class PMLDamper(object):
    """
    Class that handles the damping of the E and B fields in the PML cells,
    either on CPU or GPU.

    In agreement wih the basic principles of the PML, the damping is
    anistropic: only certain components of the fields are damped.
    """

    def __init__( self, n_pml, cdt_over_dr ):
        """
        Initialize a damping object.

        Parameters
        ----------
        # TODO
        """
        # Register the number of PML cells
        self.n_pml = n_pml

        # Create the damping arrays
        self.damp_array = generate_pml_damp_array( n_pml, cdt_over_dr )

        # Transfer the damping array to the GPU
        if cuda_installed:
            self.d_damp_array = cuda.to_device( self.damp_array )


    def damp_pml_EB( self, interp ):
        """
        Damp the fields E and B in the PML cells.

        Parameters
        ----------
        interp: list of InterpolationGrid objects (one per azimuthal mode)
            Objects that contain the fields to be damped.
        """
        # Damp the fields on the CPU or the GPU
        if interp[0].use_cuda:
            # Damp the fields on the GPU
            dim_grid, dim_block = cuda_tpb_bpg_2d( interp[0].Nz, self.n_pml )
            for m in range(len(interp)):
                cuda_damp_pml_EB[dim_grid, dim_block](
                    interp[m].Et, interp[m].Et_pml, interp[m].Ez,
                    interp[m].Bt, interp[m].Bt_pml, interp[m].Bz,
                    self.d_damp_array, self.n_pml )

        else:
            # Damp the fields on the CPU
            n_pml = self.n_pml
            for m in range(len(interp)):
                # Substract the theta PML fields to the regular theta fields
                interp[m].Et[:,-n_pml:] -= interp[m].Et_pml[:,-n_pml:]
                interp[m].Bt[:,-n_pml:] -= interp[m].Bt_pml[:,-n_pml:]
                # Damp the theta PML fields
                interp[m].Et_pml[:,-n_pml:] *= self.damp_array[np.newaxis, :]
                interp[m].Bt_pml[:,-n_pml:] *= self.damp_array[np.newaxis, :]
                # Add the theta PML fields back to the regular theta fields
                interp[m].Et[:,-n_pml:] += interp[m].Et_pml[:,-n_pml:]
                interp[m].Bt[:,-n_pml:] += interp[m].Bt_pml[:,-n_pml:]
                # Damp the z fields
                interp[m].Bz[:,-n_pml:] *= self.damp_array[np.newaxis, :]
                interp[m].Ez[:,-n_pml:] *= self.damp_array[np.newaxis, :]


def generate_pml_damp_array( n_pml, cdt_over_dr ):
    """
    Create a 1d damping array of length n_pml.

    Parameters
    ----------
    # TODO

    Returns
    -------
    A 1darray of doubles, of length n_pml, which represents the damping.
    """
    # Normalized coordinate along the pml
    x_pml = np.arange( n_pml ) * 1./n_pml

    # Calculate the damping coefficients: the damping depends on dt and dr
    # TODO: add more explanations
    damping_array = np.exp( -4. * cdt_over_dr * x_pml**2 )

    return( damping_array )


if cuda_installed:
    @cuda.jit
    def cuda_damp_pml_EB( Et, Et_pml, Ez, Bt, Bt_pml, Bz,
                      damp_array, n_pml ) :
        """
        Damp the E and B fields in the PML cells (i.e. the last n_pml cells
        in r), in an anisotropic manner which is given by the PML principles

        Parameters :
        ------------
        Et, Et_pml, Ez, Bt, Bt_pml, Bz : 2darrays of complexs
            Contain the fields to be damped
            The first axis corresponds to z and the second to r

        damp_array: 1darray of floats
            An array of length n_guards, which contains the damping factors

        n_pml: int
            Number of PML cells
        """
        # Obtain Cuda grid
        iz, i_pml = cuda.grid(2)

        # Obtain the size of the array along z and r
        Nz, Nr = Et.shape

        # Modify the fields
        if i_pml < n_pml:
            # Apply the damping arrays
            if iz < Nz:
                # Get the damping factor
                damp_factor= damp_array[i_pml]
                # Get the index in the bigger field array
                ir = Nr - n_pml + i_pml
                # Substract the theta PML fields to the regular theta fields
                Et[iz,ir] -= Et_pml[iz,ir]
                Bt[iz,ir] -= Bt_pml[iz,ir]
                # Damp the theta PML fields
                Et_pml[iz,ir] *= damp_factor
                Bt_pml[iz,ir] *= damp_factor
                # Add the theta PML fields back to the regular theta fields
                Et[iz,ir] += Et_pml[iz,ir]
                Bt[iz,ir] += Bt_pml[iz,ir]
                # Damp the z fields
                Ez[iz,ir] *= damp_factor
                Bz[iz,ir] *= damp_factor