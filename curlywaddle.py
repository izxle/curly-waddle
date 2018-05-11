#!/bin/env python

import logging

from os import path, mkdir, chdir

from ase.atoms import Atoms
from ase.io import read

from auto import gen_structs
from vasp import Vasp
from vasp.exceptions import VaspSubmitted, VaspQueued

from reader import get_args, read_config

# from typing import

# ch = logging.StreamHandler()
# ch.setLevel(logging.ERROR)
# ch.setFormatter(formatter)
# logging.basicConfig(filemode='w', style='{')
logger = logging.getLogger('curlywaddle')
# log_info = logging.getLogger('log_info')
# log_debug = logging.getLogger('log_debug')
# logger.addHandler(ch)

write_args_vasp = dict(format='vasp',
                       direct=True,
                       vasp5=True,
                       sort=True)


def run_all(slab: Atoms, adsorbate: Atoms, config: dict, **opts):
    """
        Write input files to run all combinations of slab + adsorbate
        :param slab: Atoms
        :param adsorbate: Atoms
        :param config: Dict,
        :param opts: gen_struct options
        :return:
        """
    results = []

    for ID, struct in gen_structs(slab, adsorbate, **opts):

        calc = Vasp(ID, atoms=struct, **config)
        logger.info(f'calculator for {ID} created')
        try:
            logger.info(f'{ID}: getting potential energy')
            toten = calc.potential_energy
            results.append(toten)
            # TODO: do something

        except (VaspSubmitted, VaspQueued) as e:
            logger.info(f"Couldn't get energy:\n{e}")
            results.append(None)
            print(e)

    all_jobs = len(results)
    jobs_finished = all_jobs - results.count(None)
    state = f'{jobs_finished}/{all_jobs} jobs finished'
    return state


def write_all(slab: Atoms, adsorbate: Atoms, **opts):
    """
    Write POSCAR of all combinations of slab + adsorbate
    :param slab: Atoms
    :param adsorbate: Atoms
    :param opts: gen_struct options
    :return:
    """
    results = []
    for ID, struct in gen_structs(slab, adsorbate, **opts):
        # create directory if it doesn't exist
        if not path.isdir(ID):
            mkdir(ID)
            logger.debug(f'directory {ID} created')

        fname = path.join(ID, 'POSCAR.vasp')
        try:
            struct.write(fname, **write_args_vasp)
            logger.info(f'{ID}/POSCAR written')
            results.append(1)

        except Exception as e:
            logger.error(e)
            results.append(None)

    n_structs = len(results)
    n_written_structs = n_structs - results.count(None)
    state = f'{n_written_structs}/{n_structs} structures written'
    return state


def get_adsorbate(adsorbate: str) -> Atoms:
    if adsorbate == 'O':
        ads = Atoms('O')
    elif adsorbate == 'OH':
        ads = Atoms('OH', positions=[(0, 0, 0), (0, 0, 0.96)])
    else:
        raise ValueError(f'adsorbate "{adsorbate}" not recognized.')
    return ads


def main(argv=''):
    args = get_args(argv)

    if not args.poscar_only:
        config = read_config(args.config)
    else:
        config = dict()

    # TODO: add options for gen_midpoints
    opts = dict()

    slab = read(args.slab, format=args.format)
    logger.info(f'{args.slab} Atoms object created')

    adsorbate = get_adsorbate(args.ads)
    logger.info(f'{args.ads} Atoms object created')

    if args.poscar_only:
        state = write_all(slab, adsorbate, **opts)
    else:
        state = run_all(slab, adsorbate, config, **opts)
    print(state)


if __name__ == '__main__':
    directory = path.expanduser("~/OneDrive/Documents/TAMU/structs/PtCu/111/test")
    chdir(directory)
    main('PtCu_111.vasp -p --debug')

