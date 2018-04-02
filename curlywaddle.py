#!/bin/env python

import re
from argparse import ArgumentParser
from configparser import ConfigParser
from os import path

from ase.atoms import Atoms
from ase.io import read
from vasp import Vasp
from vasp.exceptions import VaspSubmitted, VaspQueued

from auto import gen_structs
# from typing import


def get_adsorbate(adsorbate: str) -> Atoms:
    if adsorbate == 'O':
        ads = Atoms('O')
    elif adsorbate == 'OH':
        ads = Atoms('OH', positions=[(0, 0, 0), (0, 0, 0.96)])
    else:
        raise ValueError(f'adsorbate "{adsorbate}" not recognized.')
    return ads


def parse_config(config):
    debug = config.getboolean('GENERAL', 'debug')
    xc = config.get('GENERAL', 'xc')
    kpts_text = config.get('GENERAL', 'kpoints')
    kpoints = [map(int, point.split())
               for point in kpts_text.strip().split('\n')]

    incar = dict()
    for name, value in config['INCAR'].items():
        vals = []
        for v in re.split('\s*,?\s+', value):
            if re.match('-?\d+$', v):
                v = int(v)
            elif re.match('(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$', v):
                v = float(v)
            # vasp incar uses `.` around booleans, e.g. .False.
            elif v.lower().strip('. ') in ['yes', 'no', 'on', 'off', 'true', 'false', '1', '0']:
                v = v.lower().strip('. ') in ('yes', 'on', 'true', '1')
            vals.append(v)
        if len(vals) == 1:
            vals = vals[0]
        incar[name] = vals

    res = dict(kpoints=kpoints,
               xc=xc,
               debug=debug,
               **incar)
    return res


def read_config(fname: str):
    defaults = dict(GENERAL={'debug': False,
                             'xc': 'pbe'})

    parser = ConfigParser(allow_no_value=True)
    parser.read_dict(defaults)
    parser.read(fname)
    config = parse_config(parser)
    return config


def get_args(argv=''):
    parser = ArgumentParser()
    parser.add_argument('slab')
    parser.add_argument('ads', nargs='?', default='O', choices=['O', 'OH'])
    parser.add_argument('-t', '--template',
                        help='directory with template files')
    parser.add_argument('-f', '--format')
    parser.add_argument('-c', '--config', default='config.ini')
    if argv:
        if isinstance(argv, str):
            argv = argv.split()
        elif not hasattr(argv, '__iter__'):
            raise TypeError(f'argv must be `str` or iterable, not {type(argv)}')
        args = parser.parse_args(argv)
    else:
        # get arguments from terminal
        args = parser.parse_args()

    assert path.isfile(args.config), f'config: `{args.config}` must be an existing file'

    return args


def main(argv=''):
    args = get_args(argv)

    config = read_config(args.config)
    # TODO: add options for gen_midpoints
    opts = dict()

    slab = read(args.slab, format=args.format)
    adsorbate = get_adsorbate(args.ads)

    results = []

    for ID, struct in gen_structs(slab, adsorbate, **opts):

        calc = Vasp(ID, atoms=struct, **config)

        try:
            toten = calc.potential_energy
            results.append(toten)
            # TODO: do something

        except (VaspSubmitted, VaspQueued) as e:
            results.append(None)
            print(e)

    print(f'{sum(results)}/{len(results)} jobs finished')

if __name__ == '__main__':
    main()
