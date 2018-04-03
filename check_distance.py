#!/bin/env python

from argparse import ArgumentParser
from os import path, listdir

import numpy as np
from ase.geometry import distance
from ase.io import read


def atom_distance(a0, a1):
    xy0 = a0.position[:2]
    xy1 = a1.position[:2]
    return np.linalg.norm(xy1 - xy0)


def get_args(argv=''):
    parser = ArgumentParser()
    parser.add_argument('-a', '--all', action='store_true')
    parser.add_argument('--xy',
                        help='distance in the x and y axis of the given element')
    if argv:
        if isinstance(argv, str):
            argv = argv.split()
        elif not hasattr(argv, '__iter__'):
            raise TypeError(f'argv must be `str` or iterable, not {type(argv)}')
        args = parser.parse_args(argv)
    else:
        # get arguments from terminal
        args = parser.parse_args()

    return args


def main(argv='', threshold=0.001):
    args = get_args(argv)

    distances = dict()
    print('incorrect structures:')
    for name in listdir('.'):
        if not path.isdir(name):
            continue
        # get atoms object from POSCAR
        poscar_path = path.join(name, 'POSCAR')
        poscar = read(poscar_path)
        # get atoms object from CONTCAR
        contcar_path = path.join(name, 'CONTCAR')
        contcar = read(contcar_path)

        if args.xy:
            # get the first atom whose symbol == args.xy
            a0 = next(a for a in poscar if a.symbol == args.xy)
            a1 = next(a for a in contcar if a.symbol == args.xy)
            d = atom_distance(a0, a1)
        else:
            d = distance(poscar, contcar)
        #
        if d > threshold:
            print(f'{name+":":11} {d:11.3f}')
        distances[name] = d

    if args.all:
        print('all distances:')
        for name, d in distances.items():
            print(f'{name+":":11} {d:11.6f}')


if __name__ == '__main__':
    main()
