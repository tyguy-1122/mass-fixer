import json
from source.utils import calc_sequence_mass
from source.delta_finder import find_possible_deltas
# from source import delta_finder

def test_find_possible_deltas_deletions_truncations_only():
    # Load needed data
    residue_masses = json.loads(open('data/aa_masses.json').read())
    termini_species_masses = json.loads(open('data/termini_species_masses.json').read())

    # Case 1 - Easy, one possible
    desired_sequence = 'AAAA'
    observed_sequence = 'AAA'

    desired_mass = calc_sequence_mass(desired_sequence, 'H', 'OH', residue_masses, termini_species_masses)
    observed_mass = calc_sequence_mass(observed_sequence, 'H', 'OH', residue_masses, termini_species_masses)

    delta = desired_mass - observed_mass

    deltas = find_possible_deltas(desired_sequence, delta, residue_masses, 1)

    # print(deltas)

    # Case 2 - Medium 
    
def test_find_possible_deltas_long_sequence():
    print('here')