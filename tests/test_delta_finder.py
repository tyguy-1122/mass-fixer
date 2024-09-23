import json
from source.delta_finder import get_solutions
import signal
from source.peptide import Peptide
from source.delta import Delta, DeltaType
from source.delta_set import DeltaSet

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException("Method took too long to run")

def test_get_solutions_easy_deletions_truncations_only():
    # Case 1 - one possible 
    desired_sequence = 'ACDEFGHIKLMNPQRSTVWY'
    observed_sequence = 'ACDEFGHIKLNPQRSTVWY' # Missing Met

    desired_peptide = Peptide(desired_sequence, 'H', 'OH')
    observed_peptide = Peptide(observed_sequence, 'H', 'OH')

    target_mass = desired_peptide.mass - observed_peptide.mass

    solutions = get_solutions(desired_peptide, target_mass, 1)

    assert len(solutions) == 1
    assert abs(solutions[0].mass - target_mass) <= .01
    assert 'M' in [delta.description for delta in solutions[0].deltas]


    # Case 2 - three possible 
    desired_sequence = 'ACDEFGHIKLMNPQRSTVWY'
    observed_sequence = 'ACDEFGHKLMNPQRSTVWY' # Missing Ile (same mass as Leu, close to Asn)

    desired_peptide = Peptide(desired_sequence, 'H', 'OH')
    observed_peptide = Peptide(observed_sequence, 'H', 'OH')

    target_mass = desired_peptide.mass - observed_peptide.mass

    solutions = get_solutions(desired_peptide, target_mass, 1)

    assert len(solutions) == 3
    for solution in solutions:
        assert abs(solution.mass - target_mass) <= 1.01 # Greater variablity because Asn within confidence interval of Ile

    # Case 2 - one possible, multiple residues
    desired_sequence = 'ACDEFGHIKLMNPQRSTVWY'
    observed_sequence = 'ACDEFGHIKLMNQRSTWY' # Missing Pro, Val

    desired_peptide = Peptide(desired_sequence, 'H', 'OH')
    observed_peptide = Peptide(observed_sequence, 'H', 'OH')

    target_mass = desired_peptide.mass - observed_peptide.mass

    solutions = get_solutions(desired_peptide, target_mass, 1)

    assert len(solutions) == 1
    assert abs(solutions[0].mass - target_mass) <= .01
    assert 'P' in [delta.description for delta in solutions[0].deltas]
    assert 'V' in [delta.description for delta in solutions[0].deltas]

def test_find_possible_deltas_long_sequence():
    # Case 1 - 500 residues
    desired_sequence = open('resources/long-sequence-expected.txt').read()
    observed_sequence = open('resources/long-sequence-observed.txt').read() # Missing 5 random residues

    desired_peptide = Peptide(desired_sequence, 'H', 'OH')
    observed_peptide = Peptide(observed_sequence, 'H', 'OH')

    target_mass = desired_peptide.mass - observed_peptide.mass

    try: 
        signal.signal(signal.SIGALRM, timeout_handler)
        # signal.alarm(30)

        # deltas = get_solutions(desired_peptide, target_mass, 1)

        signal.alarm(0) # Disable alarm if function finishes on time

    except TimeoutException:
        assert 1 == 0
