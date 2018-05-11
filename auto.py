#!/bin/env python

from itertools import combinations
from typing import Tuple, Iterable, Iterator

import numpy as np
from ase.atoms import Atoms
from ase.geometry import get_layers
from numpy.linalg import  norm

import logging

logger = logging.getLogger('curlywaddly')


class MyAtoms(Atoms):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.midpoints = []

    def in_cell(self, index):
        return all(0 <= coord < 1
                   for coord in self.get_scaled_positions()[index])

    def calc_midpoints(self):
        # TODO: calc midpoints
        midpoints = []
        self.midpoints = midpoints

    def get_midpoints(self):
        if not self.midpoints:
            self.calc_midpoints()
        return self.midpoints


def in_cell(atoms: Atoms, index=None, position=None):
    '''
    Checks if `x` and `y` scaled coordinates are between 0 and 1
    :param atoms: MUST have scaled_positions property
    :param index: atom index
    :param position:
    :return: bool
    '''
    if index is not None:
        scaled = atoms.get_scaled_positions()[index][:2]
    elif position is not None:
        # scales position to cell
        scaled = np.linalg.solve(atoms.cell.T, position.T).T[:2]
    else:
        raise ValueError('at least one of `position` or `index` '
                         'must be provided')
    return all(0 <= coord < 1
               for coord in scaled)


def valid_dist(indices, atoms, cutoff_dist):
    """
    Return whether
    :param indices:
    :param atoms:
    :param cutoff_dist:
    :return:
    """
    return all(atoms.get_distance(i, j) < cutoff_dist
               for i, j in combinations(indices, 2))


def valid_comb(indices, atoms, expanded_atoms, cutoff_dist):
    return any(in_cell(atoms, position=expanded_atoms[index].position) for index in indices) \
           and valid_dist(indices, expanded_atoms, cutoff_dist)


def distance(point1, point2):
    return norm(point1 - point2)


def get_incenter(point1, point2, point3):
    Ax, Ay, Az = point1
    Bx, By, Bz = point2
    Cx, Cy, Cz = point3
    a = distance(point2, point3)
    b = distance(point1, point3)
    c = distance(point1, point2)
    p = a + b + c
    Ox = (a*Ax + b*Bx + c*Cx) / p
    Oy = (a*Ay + b*By + c*Cy) / p
    Oz = np.mean( (Az, Bz, Cz) )
    return np.array( (Ox, Oy, Oz) )


def gen_midpoints(slab: Atoms, cutoff_dist: float=3.5, heights: Iterable[float]=None, **kw):
    """

    :param slab:
    :param cutoff_dist:
    :param heights:
    :param kw:
    :return:
    """
    # default heights
    h = heights or [(0, 0, i) for i in (0, 2, 1.8, 1.5, 1.3)]
    tags, layer_pos = get_layers(slab, (0, 0, 1), 0.3)
    surface_atoms = slab[tags == max(tags)]
    # n_atoms = len(surface_atoms)
    atoms = surface_atoms.repeat([3, 3, 1])
    n_points = len(atoms)
    # TODO: check if atoms moved
    atoms.wrap(center=(1/6, 1/6, 0))

    mem = set()

    def memo(pos: Tuple[float, float]) -> bool:
        label = f'{pos[0]:.2f} {pos[1]:.2f}'
        res = label in mem
        if not res:
            mem.add(label)
        return res

    for i, pos in enumerate(surface_atoms.positions):
        memo(pos)
        yield (f'1_{i}', pos + h[1])

    for i in range(2, len(h)):
        ix = 0
        for indices in combinations(range(n_points), i):
            if not valid_comb(indices, surface_atoms, atoms, cutoff_dist):
                continue
            positions = atoms.positions[list(indices)]
            # TODO:
            if i == 2:
                mean = np.mean(positions, axis=0)
            else:
                mean = get_incenter(positions[0], positions[1], positions[2])
            mean += h[i]
            if not in_cell(surface_atoms, position=mean) or memo(tuple(mean)):
                continue
            yield (f'{i}_{ix}', mean)
            ix += 1


def center_at_origin(atoms: Atoms) -> Atoms:
    """
    Returns a copy of the atoms object centered at the lowest atoms in the z axis
    :param atoms: Atoms object
    :return: Copy of the translated Atoms object
    """
    ads = atoms.copy()
    if len(ads) == 1:
        ads[0].position = (0, 0, 0)
    else:
        tags, layer_pos = get_layers(ads, (0, 0, 1), 0.3)

        if sum(tags == 0) == 1:
            center = atoms[0].position
        else:
            lowest_atoms = ads[tags == 0]
            center = np.mean([a.position for a in lowest_atoms], axis=0)

        ads.translate(-center)
    return ads


def gen_structs(atoms: Atoms, adsorbate: Atoms, **kw) -> Iterator[Tuple[str, Atoms]]:
    """

    :param atoms:
    :param adsorbate:
    :param kw:
    :return: Tuple[ID, atoms + adsorbate]
    """
    # TODO: orient/rotate
    ads = center_at_origin(adsorbate)

    for ID, point in gen_midpoints(atoms, **kw):
        a = ads.copy()
        a.translate(point)
        a.set_scaled_positions(point)
        slab_ads = atoms + a
        yield (ID,  slab_ads)

