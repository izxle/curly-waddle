#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from itertools import count
from ase.io import read
from ase.atoms import Atoms
from ase.calculators.vasp import Vasp
from typing import Iterable
from numbers import Real
import logging

# single point calculations of small systems should be fast


def run_vasp(atoms: Atoms):
    xc = 'PBE'
    init_kpoints = (5, 5, 5)
    incar = dict(
        isif=3,
        ibrion=1,
        sigma=0.05,
        nelm=200,
        nsw=200)
    calc = Vasp(xc=xc,
                kpoints=init_kpoints,
                **incar)
    atoms.set_calculator(calc)

    plot_params = dict(save_plot=True)
    encut_params = dict(start=250, step=50)
    optimize(atoms, 'encut', params=encut_params,
             threshold=1e-2, plot_params=plot_params)

    # kpoint_params = dict(start=3, step=1)
    # optimize(atoms, 'kpoints', params=kpoint_params,
    #          threshold=1e-2, plot_params=plot_params)
    return


def gen_calculation(atoms: Atoms, parameter: str, start: int, step: int=1):
    for val in count(start, step):
        atoms.calc.set(**{parameter: val})
        yield (val, atoms.get_total_energy())


def optimize(atoms: Atoms, parameter: str, params: dict={},
             threshold: float=1e-2, max_iterations: int=15,
             plot_params: dict={}):
    values = []
    energies = []
    # TODO: add optimization algorithm
    # naive approach
    prev_energy = float('inf')
    for i, (value, energy) in enumerate(gen_calculation(atoms, parameter, **params)):
        values.append(value)
        energies.append(energy)
        if abs(energy - prev_energy) < threshold or i >= max_iterations:
            break
        prev_energy = energy

    # plotting section
    if plot_params:
        plot(parameter, values, energies, **plot_params)


def plot(name: str, x: Iterable[Real], y: Iterable[Real],
         ax: Axes=None, save_plot: bool=False, **kwargs):
    if ax is None:
        fig = plt.figure(name)
        ax = fig.add_subplot(111)

    ax.plot(x, y, **kwargs)

    if save_plot:
        plt.savefig(f'{name}.png')


def get_args(argv: str=''):
    """
    argument parser for command line execution
    :param argv: string emulating the command line arguments
    :return: arguments Namespace
    """
    # initialization
    kw = dict(description='',
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser = argparse.ArgumentParser(**kw)

    # argument definition
    parser.add_argument('atoms',
                        help='filename to read atoms from')
    parser.add_argument('-f', '--format',
                        help='format of atoms file')

    # read arguments
    if isinstance(argv, str):
        argv = argv.split()
    elif not hasattr(argv, '__iter__'):
        raise TypeError(f'argv must be `str` or iterable, not {type(argv)}')
    args = parser.parse_args(argv) if argv else parser.parse_args()

    return args


def main(argv=''):
    args = get_args(argv)
    atoms = read(args.atoms, format=args.format)
    run_vasp(atoms)


if __name__ == '__main__':
    main()
