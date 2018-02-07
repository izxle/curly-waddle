#!/bin/env python

from itertools import combinations
from typing import Tuple, Iterable, Iterator

import numpy as np
from ase.atoms import Atoms
from ase.geometry import get_layers


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


def valid_comb(points, atoms, cutoff_dist):
    return any(in_cell(atoms, index=point) for point in points) \
           and valid_dist(points, atoms, cutoff_dist)


def gen_midpoints(slab: Atoms, cutoff_dist: float=3.5, heights: Iterable[float]=None, **kw):
    h = heights or [(0, 0, i) for i in (0, 2, 1.8, 1.5, 1.3)]
    tags, layer_pos = get_layers(slab, (0, 0, 1), 0.3)
    surface_atoms = slab[tags == max(tags)]
    # n_atoms = len(surface_atoms)
    atoms = surface_atoms.repeat([3, 3, 1])
    n_points = len(atoms)
    # TODO: check if atoms moved
    atoms.wrap(center=(0, 0, 0))

    mem = set()

    def memo(pos):
        label = f'{pos[0]:.2f} {pos[1]:.2f}'
        res = label in mem
        if not res:
            mem.add(label)
        return res

    # combs = {f'1_{i}': pos + h[1] for i, pos in enumerate(surface_atoms.positions)}
    for i, pos in enumerate(surface_atoms.positions):
        memo(pos)
        yield (f'1_{i}', pos + h[1])
    for i in range(2, len(h)):
        ix = 0
        for comb in combinations(range(n_points), i):
            if not valid_comb(comb, atoms, cutoff_dist):
                continue
            mean = np.mean(comb, axis=0) + h[i]
            if not in_cell(atoms, position=mean) or memo(mean):
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
    # yield from ((ID, atoms + ads.copy().translate(point))
    #             for ID, point in gen_midpoints(atoms, **kw))
    for ID, point in gen_midpoints(atoms, **kw):
        a = ads.copy()
        a.translate(point)
        yield (ID, atoms + a)


# def getArgs():
#     parser = ArgumentParser()
#     parser.add_argument('filename')
#     parser.add_argument('adsorbate', dest='ads', nargs='?', default='O',
#                         choices=['O', 'OH'])
#     parser.add_argument('-f', '--format')
#     parser.add_argument('--template',
#                         help='folder with template files')
#     return parser.parse_args()


# def get_adsorbate(adsorbate):
#     if adsorbate == 'O':
#         ads = Atoms('O')
#     elif adsorbate == 'OH':
#         ads = Atoms('OH', positions=[(0, 0, 0), (0, 0, 0.96)])
#     return ads

# def main():
#     args = getArgs()
#     templates = listdir(args.template)
#     atoms = read(args.filename, args.format) if args.format else \
#         read(args.filename)
#     ads = get_adsorbate(args.ads)
#
#     for ID, struct in gen_structs(atoms, ads, **vars(args)):
#         exists(ID) or mkdir(ID)
#         for tmp in templates:
#             copy(join(args.template, tmp), join(ID, tmp))
#         struct.write(join(ID, 'POSCAR'), format='vasp', vasp5=True,
#                      sort=True, direct=True)
