#!/bin/env python

from argparse import ArgumentParser
from os import path, listdir
from ase.io import read
from ase.geometry import distance

def get_args(argv=''):
    parser = ArgumentParser()
    parser.add_argument('something')
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

    distances = list()
    for name in listdir('.'):
        if not path.isdir(name):
            continue
        rel_path_0 = path.join(name, 'POSCAR')
        rel_path_1 = path.join(name, 'CONTCAR')

        a0 = read(rel_path_0)
        a1 = read(rel_path_1)

        d = distance(a0, a1)

        distances.append(d)

    for d in distances:
        if d < threshold:
            pass




if __name__ == '__main__':
    main()
