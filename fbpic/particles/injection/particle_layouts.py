# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)
It defines a set of numerical particle layout.
"""
import warnings

class PseudoRandomLayout( object ):

    def __init__(self, n_macroparticles=None, seed=None, **kw):
        """
        TODO
        """
        assert n_macroparticles is not None
        self.n_macroparticles = n_macroparticles
        if seed is not None:
            warnings.warn('Fixing the seed is not supported in fbpic.')

class GriddedLayout( object ):
    """
    TODO
    """
    def __init__(self, n_macroparticle_per_cell, **kw):
        self.n_macroparticle_per_cell = n_macroparticle_per_cell
        assert 'r' in n_macroparticle_per_cell.keys()
        assert 'z' in n_macroparticle_per_cell.keys()
        assert 'theta' in n_macroparticle_per_cell.keys()
