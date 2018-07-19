# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)
It defines a set of numerical particle layout.
"""
from lpa_utils.boost import BoostConverter

# Particle layout classes
# -----------------------

class PseudoRandomLayout( object ):

    def __init__(self, n_macroparticles=None, seed=None, **kw):
        """
        TODO
        """
        assert n_macroparticles is not None
        self.n_macroparticles = n_macroparticles
        assert seed is None # Seed is currently not supported

    def generate_particles( self, distribution, comm, gamma_boost ):
        """
        TODO
        """
        if not hasattr( distribution, generate_random_particles ):
            raise ValueError('')
        x, y, z, ux, uy, uz, w = \
            distribution.generate_random_particles( self.n_macroparticles )

        x, y, z, ux, uy, uz, inv_gamma, w = \
            boost_and_select( x, y, z, ux, uy, uz, w, comm, gamma_boost )

        return len(w), x, y, z, ux, uy, uz, inv_gamma, w
