#!/bin/env python

from argparse import ArgumentParser, Namespace
from os import mkdir, listdir
from os.path import exists, join
from ase.atoms import Atoms
from ase.io import read
# from typing import

from auto import gen_structs


def get_adsorbate(adsorbate: str) -> Atoms:
    if adsorbate == 'O':
        ads = Atoms('O')
    elif adsorbate == 'OH':
        ads = Atoms('OH', positions=[(0, 0, 0), (0, 0, 0.96)])
    else:
        raise ValueError(f'adsorbate "{adsorbate}" not recognized.')
    return ads


def get_args(argv=''):
    parser = ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('ads', nargs='?', default='O', choices=['O', 'OH'])
    parser.add_argument('-t', '--template',
                        help='directory with template files')
    parser.add_argument('-f', '--format')
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


def main(argv=''):
    args = get_args(argv)

    atoms = read(args.filename, format=args.format)
    adsorbate = get_adsorbate(args.ads)

    for ID, struct in gen_structs(atoms, adsorbate, **vars(args)):
        if not exists(ID):
            mkdir(ID)
        # for tmp in templates:
        #     copy(join(args.template, tmp), join(ID, tmp))
        struct.write(join(ID, 'POSCAR'), format='vasp', vasp5=True,
                     sort=True, direct=True)


if __name__ == '__main__':
    main('../VaspTools/POSCAR.draft -f vasp')
