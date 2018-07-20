# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)
It defines a set of numerical particle layout.
"""

class PseudoRandomLayout( object ):

    def __init__(self, n_macroparticles=None, seed=None, **kw):
        """
        TODO
        """
        assert n_macroparticles is not None
        self.n_macroparticles = n_macroparticles
        assert seed is None # Seed is currently not supported
