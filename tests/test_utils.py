import json
from source.utils import calc_sequence_mass

def test_calc_sequence_mass_full_length_natural_peptide():
    aa_masses = json.loads(open('data/aa_masses.json').read())
    termini_masses = json.loads(open('data/termini_species_masses.json').read())

    # Short peptide
    true_mass = 380.44
    calculated_mass = calc_sequence_mass('KAY', 'H', 'OH', aa_masses, termini_masses)

    assert abs(true_mass - calculated_mass) <= .01

    # Long
    true_mass = 2296.5973
    calculated_mass = calc_sequence_mass('ACDEFGHIKLMNPQRSTWY', 'H', 'OH', aa_masses, termini_masses)

    assert abs(true_mass - calculated_mass) <= .01

def test_calc_sequence_mass_special_termini():
    aa_masses = json.loads(open('data/aa_masses.json').read())
    termini_masses = json.loads(open('data/termini_species_masses.json').read())

    # Fmoc, NHNH2
    true_mass = 615.75
    calculated_mass = calc_sequence_mass('KAY', 'Fmoc', 'NHNH2', aa_masses, termini_masses)

    assert abs(true_mass - calculated_mass) <= .01

    # Acetyl N-cap, MPAA
    true_mass = 571.4948
    calculated_mass = calc_sequence_mass('KAY', 'Acetyl N-cap', 'MPAA', aa_masses, termini_masses)

    assert abs(true_mass - calculated_mass) <= .01

    # H, NH2
    true_mass = 379.4578
    calculated_mass = calc_sequence_mass('KAY', 'H', 'NH2', aa_masses, termini_masses)

    assert abs(true_mass - calculated_mass) <= .01

    # H, Dbz
    true_mass = 513.5058
    calculated_mass = calc_sequence_mass('KAY', 'H', 'Dbz', aa_masses, termini_masses)

    assert abs(true_mass - calculated_mass) <= .01

    # H, Nbz
    true_mass = 539.5878
    calculated_mass = calc_sequence_mass('KAY', 'H', 'Nbz', aa_masses, termini_masses)

    assert abs(true_mass - calculated_mass) <= .01

def test_calc_sequence_mass_with_non_canoncials():
    aa_masses = json.loads(open('data/aa_masses.json').read())
    non_canonicals = json.loads(open('data/non-canonical_aas.json').read())
    for non_canonical in non_canonicals:
        aa_masses[non_canonicals[non_canonical]['symbol']] = non_canonicals[non_canonical]['mass']
    termini_masses = json.loads(open('data/termini_species_masses.json').read())

    # Short peptide
    true_mass = 1022.4653
    calculated_mass = calc_sequence_mass('ZXOUB', 'H', 'OH', aa_masses, termini_masses)

    assert abs(true_mass - calculated_mass) <= .01

def test_calc_sequence_mass_with_non_canoncials_and_special_termini():
    aa_masses = json.loads(open('data/aa_masses.json').read())
    non_canonicals = json.loads(open('data/non-canonical_aas.json').read())
    for non_canonical in non_canonicals:
        aa_masses[non_canonicals[non_canonical]['symbol']] = non_canonicals[non_canonical]['mass']
    termini_masses = json.loads(open('data/termini_species_masses.json').read())

    # Short peptide
    true_mass = 643.6328
    calculated_mass = calc_sequence_mass('KAYZ', 'H', 'MPAA', aa_masses, termini_masses)

    assert abs(true_mass - calculated_mass) <= .01