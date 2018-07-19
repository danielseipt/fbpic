# Copyright 2018, FBPIC contributors
# Authors: Remi Lehe, Manuel Kirchen
# License: 3-Clause-BSD-LBNL
"""
This file is part of the Fourier-Bessel Particle-In-Cell code (FB-PIC)

It defines the picmi laser interface
"""
# Laser injector objects
# ----------------------

# TODO: should we rename these objects with `Injector` at the end?
# TODO: should we add an antenna velocity?
# TODO: should we add direct injection in PICMI
from PICMI_Base import PICMI_LaserAntenna, PICMI_GaussianLaser

class LaserAntenna( PICMI_LaserAntenna ):

    def init(self):
        if self.antenna_normal != (0,0,1.):
            raise ValueError(
            'In fbpic, the normal vector should be along z. (i.e. (0,0,1.))')
        if self.antenna_position[0] != 0 or self.antenna_position[1] != 0:
            raise ValueError(
            'In fbpic, the position of the antenna should be on the z axis.')

class LaserDirect( object ):

    def __init__(self):
        pass

# Laser profile objects
# ---------------------

class GaussianLaser( PICMI_GaussianLaser ):

    def init(self):
        if self.focal_position[0] != 0 or self.focal_position[1] != 0:
            raise ValueError(
            'In fbpic, the focal position of a laser should be on the z axis.')
        if self.centroid_position[0] != 0 or self.centroid_position[1] != 0:
            raise ValueError(
            'In fbpic, the centroid position of a laser should be on the z axis.')
        if self.propagation_direction[0] != 0 or self.propagation_direction[1] != 0:
            raise ValueError(
            'In fbpic, the propagation position of a laser should be along the z axis.')
