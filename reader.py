import logging
import re
from argparse import ArgumentParser
from configparser import ConfigParser
from os import path


logger = logging.getLogger('curlywaddle')
formatter = logging.Formatter('{levelname:8}:{filename}:{funcName}:{lineno}: {message}', style='{')


def parse_config(config) -> dict:
    """
    Parses information from configuration
    :param config: ConfigParser
    :return: dict
    """
    if config.getboolean('GENERAL', 'debug'):
        debug = 10
    else:
        debug = 0

    xc = config.get('GENERAL', 'xc')

    # KPOINTS section
    kpts_text = config.get('GENERAL', 'kpoints')
    kpoints = [list(map(int, point.split()))
               for point in kpts_text.strip().split('\n')]
    if len(kpoints) == 1:
        kpoints = kpoints[0]

    # INCAR section
    incar = dict()
    for name, value in config['INCAR'].items():
        vals = []
        for v in re.split('\s*,?\s+', value):
            if re.match('[-+]?\d+$', v):
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

    res = dict(kpts=kpoints,
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
    parser.add_argument('-f', '--format')
    parser.add_argument('-p', '--poscar_only', action='store_true',
                        help='writes POSCAR only, does not run')
    parser.add_argument('-c', '--config', default='config.ini')
    parser.add_argument('--debug', action='store_true')
    if argv:
        if isinstance(argv, str):
            argv = argv.split()
        elif not hasattr(argv, '__iter__'):
            raise TypeError(f'argv must be `str` or iterable, not {type(argv)}')
        args = parser.parse_args(argv)
    else:
        # get arguments from terminal
        args = parser.parse_args()

    if args.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logger.setLevel(log_level)
    fh = logging.FileHandler('out.log')
    fh.setFormatter(formatter)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)

    logger.debug(f'args:{vars(args)}')

    if not args.poscar_only:
        assert path.isfile(args.config), f'config: `{args.config}` must be an existing file'

    return args
