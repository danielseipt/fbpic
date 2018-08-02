# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)

It defines various elements of the PICMI interface
"""
from PICMI_Base import PICMI_CylindricalGrid, PICMI_ElectromagneticSolver

class CylindricalGrid(PICMI_CylindricalGrid):

    def init(self):
        if self.bc_rmax != 'reflective':
            raise ValueError(
            "In fbpic, `bc_rmax` can only be 'reflective' - not %s." %self.bc_rmax)
        if self.bc_zmin not in ['periodic', 'open']:
            raise ValueError("In fbpic, `bc_zmin` can only be "
                "'periodic' or 'open' - not %s." %self.bc_zmin)
        if self.bc_zmax != self.bc_zmax:
            raise ValueError(
            "In fbpic, `bc_zmax` should be identical to `bc_zmin`.")
        if self.moving_window_velocity[0] !=0 or self.moving_window_velocity[1]:
            raise ValueError(
            "In fbpic, `moving_window_velocity` should be along z.")

class ElectromagneticSolver(PICMI_ElectromagneticSolver):

    def init(self, **kw ):
        pass
