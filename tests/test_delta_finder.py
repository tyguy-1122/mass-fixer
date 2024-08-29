import json
from source.utils import calc_sequence_mass
from source.delta_finder import find_possible_deltas
import signal


class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException("Method took too long to run")

def test_find_possible_deltas_easy_deletions_truncations_only():
    # Load needed data
    residue_masses = json.loads(open('data/aa_masses.json').read())
    termini_species_masses = json.loads(open('data/termini_species_masses.json').read())

    # Case 1 - one possible 
    desired_sequence = 'ACDEFGHIKLMNPQRSTVWY'
    observed_sequence = 'ACDEFGHIKLNPQRSTVWY' # Missing Met

    desired_mass = calc_sequence_mass(desired_sequence, 'H', 'OH', residue_masses, termini_species_masses)
    observed_mass = calc_sequence_mass(observed_sequence, 'H', 'OH', residue_masses, termini_species_masses)

    delta = desired_mass - observed_mass

    deltas = find_possible_deltas(desired_sequence, delta, residue_masses, 1)

    assert len(deltas) == 1
    assert abs(deltas[0]['mass'] - delta) <= .01
    assert deltas[0]['residues'] == [10]


    # Case 2 - three possible 
    desired_sequence = 'ACDEFGHIKLMNPQRSTVWY'
    observed_sequence = 'ACDEFGHKLMNPQRSTVWY' # Missing Ile (same mass as Leu, close to Asn)

    desired_mass = calc_sequence_mass(desired_sequence, 'H', 'OH', residue_masses, termini_species_masses)
    observed_mass = calc_sequence_mass(observed_sequence, 'H', 'OH', residue_masses, termini_species_masses)

    delta = desired_mass - observed_mass

    deltas = find_possible_deltas(desired_sequence, delta, residue_masses, 1)

    assert len(deltas) == 3
    for solution in deltas:
        assert abs(solution['mass'] - delta) <= 1.01 # Greater variablity because Asn within confidence interval of Ile

    # Case 2 - one possible, multiple residues
    desired_sequence = 'ACDEFGHIKLMNPQRSTVWY'
    observed_sequence = 'ACDEFGHIKLMNQRSTWY' # Missing Pro, Val

    desired_mass = calc_sequence_mass(desired_sequence, 'H', 'OH', residue_masses, termini_species_masses)
    observed_mass = calc_sequence_mass(observed_sequence, 'H', 'OH', residue_masses, termini_species_masses)

    delta = desired_mass - observed_mass

    deltas = find_possible_deltas(desired_sequence, delta, residue_masses, 1)

    assert len(deltas) == 1
    assert abs(deltas[0]['mass'] - delta) <= .01
    assert 12 in deltas[0]['residues']
    assert 17 in deltas[0]['residues']

def test_find_possible_deltas_long_sequence():
    # Load needed data
    residue_masses = json.loads(open('data/aa_masses.json').read())
    termini_species_masses = json.loads(open('data/termini_species_masses.json').read())

    # Case 1 - 500 residues
    desired_sequence = open('resources/long-sequence-expected.txt').read()
    observed_sequence = open('resources/long-sequence-observed.txt').read() # Missing 5 random residues

    desired_mass = calc_sequence_mass(desired_sequence, 'H', 'OH', residue_masses, termini_species_masses)
    observed_mass = calc_sequence_mass(observed_sequence, 'H', 'OH', residue_masses, termini_species_masses)

    delta = desired_mass - observed_mass

    try: 
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(30)

        deltas = find_possible_deltas(desired_sequence, delta, residue_masses, 1)

        signal.alarm(0) # Disable alarm if function finishes on time

    except TimeoutException:
        assert 1 == 0
