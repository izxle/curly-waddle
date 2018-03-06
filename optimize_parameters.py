#!/usr/bin/env python
import argparse
import ase


# single point calculations should be fast


def get_args(argv: str=''):
    """
    argument parser for command line execution
    :param argv: string emulating the command line arguments
    :return: arguments Namespace
    """
    kw = dict(description='',
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser = argparse.ArgumentParser(**kw)

    # argument definition
    parser.add_argument('first_argument')

    if isinstance(argv, str):
        argv = argv.split()
    elif not hasattr(argv, '__iter__'):
        raise TypeError(f'argv must be `str` or iterable, not {type(argv)}')
    args = parser.parse_args(argv) if argv else parser.parse_args()

    return args


def main(argv=''):
    args = get_args(argv)


if __name__ == '__main__':
    main()
