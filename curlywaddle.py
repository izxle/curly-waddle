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


def get_args():
    parser = ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('adsorbate', nargs='?', default='O',
                        choices=['O', 'OH'])
    parser.add_argument('-f', '--format')
    return parser.parse_args()


def main():
    args = get_args()

    atoms = read(args.filename, format=args.format)
    ads = get_adsorbate(args.adsorbate)

    for ID, struct in gen_structs(atoms, ads, **vars(args)):
        if not exists(ID):
            mkdir(ID)
        # for tmp in templates:
        #     copy(join(args.template, tmp), join(ID, tmp))
        struct.write(join(ID, 'POSCAR'), format='vasp', vasp5=True,
                     sort=True, direct=True)


if __name__ == '__main__':
    main()