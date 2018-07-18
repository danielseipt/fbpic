# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)

It defines the picmi interface
"""
# Define general variables that each PICMI code should define
code = 'fbpic'

# Import object defined in other files
from .grid import CylindricalGrid, ElectromagneticSolver, BinomialSmoother
from .simulation import Simulation

__all__ = ['code', 'CylindricalGrid', 'ElectromagneticSolver', 'Simulation', \
            'BinomialSmoother' ]
