
from ase.io import write
from ase.constraints import FixAtoms, FixScaled
from itertools import combinations, product
import numpy as np
from ase.geometry import cell_to_cellpar, get_layers
from ase.atoms import Atoms

def adsorb_on_all(atoms, ads, cutoff_dist=3.5, symmetry=None):
    '''
    :param atoms: 
    :param ads: 
    :param reduced_cell: 
    :param cutoff_dist: 
    :return: 
    '''

    h = [None,
         [0, 0, 1.3],
         [0, 0, 1.1],
         [0, 0, 0.98],
         [0, 0, 0.95]]

    assert isinstance(atoms, Atoms), '"atoms" must be Atoms instance'
    assert isinstance(ads, Atoms), '"ads" must be Atoms instance'

    init_atoms = atoms.copy()
    ads = ads.copy()

    def in_cell(position):
        fractional = np.linalg.solve(cell.T, position.T).T
        return all(0 <= c < 1 for c in fractional)

    def in_cell_scaled(pos):
        return all(0 <= c < 1 for c in pos)

    def valid_dist_ase(indices):
        return all(atoms.get_distance(i, j) < cutoff_dist
                   for i, j in combinations(indices, 2))

    def valid_comb_ase(indices):
        return any(in_cell_scaled(scaled_pos[i]) for i in indices)\
               and valid_dist_ase(indices)

    def valid_dist(points):
        return all(abs(np.linalg.norm(p1 - p2)) < cutoff_dist
                   for p1, p2 in combinations(points, 2))

    def valid_comb(points):
        return any(in_cell(point) for point in points) and valid_dist(points)

    def get_mean_points_ase(indices):
        points = [positions[i] for i in indices]
        mean = np.mean(points, axis=0) + h[len(indices)]
        mean_points = [mean]
        mean_points += [np.mean([point, mean], axis=0) for point in points]
        return mean_points

    def valid_point(point):
        return in_cell(point) and not in_mem(point)

    def filter_valid_points(points):
        for p in points:
            if valid_point(p):
                yield p

    def get_mean_points(points):
        mean = np.mean(points, axis=0) + h[len(points)]
        mean_points = [mean]
        mean_points += [np.mean([point, mean], axis=0) for point in points]
        return mean_points

    def in_mem(point):
        coords = '{:.2f} {:.2f}'.format(*point)
        res = coords in mem
        if not res:
            mem.add(coords)
        return res

    def in_mins(point):
        return '{:.2f} {:.2f}'.format(*point) in minima

    cell = init_atoms.get_cell()

    # get atoms of last layer
    tags, layer_pos = get_layers(atoms, (0, 0, 1), 0.3)
    surface = atoms[tags == max(tags)]
    atoms = surface.repeat([3, 3, 1])
    # TODO: translate in fractional coordinates
    #atoms.translate([-a, -b, 0])
    atoms.wrap(center=(0,0,0))
    positions = atoms.get_positions()
    scaled_pos = atoms.get_scaled_positions()

    n_atoms_init = len(init_atoms)
    n_atoms = len(atoms)

    # get all points where to adsorbe adsorbate
    tops = [pos + h[1] for pos in surface.positions]
    mem = {'{:.2f} {:.2f}'.format(*pos) for pos in tops}
    combs = {f'1_{i}': pos
             for i, pos in enumerate(tops)}
    count = {i: 0 for i in range(2, len(h))}
    count[1] = len(tops)
    minima = set(mem)

    for i in range(2, 5):
        #print(f'combs {i}')
        for comb in combinations(range(n_atoms), i):
        #for j, comb in enumerate(combinations(positions, i), 1):
            #print(f'\r{i}_{j}', end='')
            if not valid_comb_ase(comb): continue
            mean_points = get_mean_points_ase(comb)
            #print(f'\r{i}_{j}: {len(mean_points)} mean points')

            minimum = mean_points[0]
            if in_cell(minimum):
                minima.add('{:.2f} {:.2f}'.format(*minimum))

            for point in filter_valid_points(mean_points):
                min_label = 'm' if in_mins(point) else 't'
                ID = f'{min_label}_{len(comb)}_{count[len(comb)]}'
                #print(f'\r{ID}', end='')
                combs[ID] = point
                count[len(comb)] += 1

            #print(f'{count[len(comb)]} valid points')
    print(f'mem: <{len(mem)}>\tcombs: <{len(combs)}>\tmins: <{len(minima)}>')

    # adsorb on all points
    if len(ads) == 1:
        ads[0].position = (0,0,0)
    else:
        tags, layer_pos = get_layers(ads, (0,0,1), 0.3)
        lowest_atoms = ads[tags==0]
        center = np.mean([a.position for a in lowest_atoms], axis=0)
        ads.translate(-center)

    structs = {}
    formula = ads.get_chemical_formula()
    #test = init_atoms.copy()

    constraint = FixAtoms(indices=range(n_atoms_init))
    constraints = [FixScaled(cell, i, (1, 1, 0))
                   for i in range(n_atoms_init, n_atoms_init + len(ads))]

    for ID, pos in combs.items():
        #print(f'{ID}')
        a = ads.copy()
        a.translate(pos)
        filename = f'POSCAR_{formula}_{ID}.vasp'
        tmp_struct = init_atoms + a
        tmp_struct.set_constraint([constraint] + constraints)
        structs[filename] = tmp_struct
        #test.extend(a)

    #with open('test.out', 'w') as f:
    #    for ID in sorted(combs.keys()):
    #        f.write(ID + '\n')

    #write_vasp('test_all.vasp', test)
    #from ase.io import write
    #write('test_atoms.cif', atoms)

    for filename, atoms in structs.items():
        write_vasp(filename, atoms)

def write_vasp(filename, atoms):
    def kw():
        format = 'vasp'
        direct = True
        vasp5 = True
        sort = True
        return locals()

    #from vasp import write_vasp as write
    write(filename, atoms, **kw())





