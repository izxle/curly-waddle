
import logging

from ase.atoms import Atoms

from vasp import Vasp
from vasp.exceptions import VaspQueued, VaspSubmitted

logger = logging.getLogger('curlywaddly')


def relax_struct(atoms: Atoms, config: dict):
    default_config=dict(encut=400,
                        isif=3,
                        ibrion=2,
                        ispin=2,
                        lreal='Auto')
    default_config.update(config)
    calc = Vasp('init_relax', atoms=atoms, **config)
    try:
        energy = calc.potential_energy
        print(energy)
        state = 'initial relaxation complete'
    except (VaspSubmitted, VaspQueued) as e:
        logger.info(f"Couldn't get energy:\n{e}")
        print(e)
        state = 'initial relaxation running'

    return state
