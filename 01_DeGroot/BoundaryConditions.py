from enum import Enum

class BoundaryLocation(Enum):
    WEST = 1    ;     EAST = 2

class DirichletBc:
    """Class defining a Dirichlet boundary condition"""

    def __init__(self, phi, grid, value, loc):
        self._phi = phi     ; self._grid = grid     ; self._value = value   ; self._loc = loc

    def value(self):                return self._value
    def coeff(self):                return 0

    def apply(self):
        if self._loc is BoundaryLocation.WEST:      self._phi[0] = self._value
        elif self._loc is BoundaryLocation.EAST:    self._phi[-1] = self._value
        else:   raise ValueError("Unknown boundary location")
        

class NeumannBc:

    def __init__(self, phi, grid, gradient, loc):
        self._phi = phi ; self._grid = grid ; self._gradient = gradient ; self._loc = loc

    def value(self):
        if self._loc is BoundaryLocation.WEST:          return self._phi[1] - self._gradient*self._grid.dx_WP[0]
        elif self._loc is BoundaryLocation.EAST:        return self._phi[-2] + self._gradient*self._grid.dx_PE[-1]
        else:                                           raise ValueError("Unknown boundary location")

    def coeff(self):
        return 1

    def apply(self):
        if self._loc is BoundaryLocation.WEST:      self._phi[0] = self._phi[1] - self._gradient*self._grid.dx_WP[0]
        elif self._loc is BoundaryLocation.EAST:    self._phi[-1] = self._phi[-2] + self._gradient*self._grid.dx_PE[-1]
        else:                                       raise ValueError("Unknown boundary location")
