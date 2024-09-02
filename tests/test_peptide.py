import json
from source.peptide import Peptide

def test_calc_sequence_mass_full_length_natural_peptide():
    # Short peptide
    true_mass = 380.44
    peptide = Peptide('KAY', 'H', 'OH')

    assert abs(true_mass - peptide.mass) <= .01

    # Long
    true_mass = 2296.5973
    peptide = Peptide('ACDEFGHIKLMNPQRSTWY', 'H', 'OH')

    assert abs(true_mass - peptide.mass) <= .01

def test_calc_sequence_mass_special_termini():
    # Fmoc, NHNH2
    true_mass = 615.75
    peptide = Peptide('KAY', 'Fmoc', 'NHNH2')

    assert abs(true_mass - peptide.mass) <= .01

    # Acetyl N-cap, MPAA
    true_mass = 571.4948
    peptide = Peptide('KAY', 'Acetyl N-cap', 'MPAA')

    assert abs(true_mass - peptide.mass) <= .01

    # H, NH2
    true_mass = 379.4578
    peptide = Peptide('KAY', 'H', 'NH2')

    assert abs(true_mass - peptide.mass) <= .01

    # H, Dbz
    true_mass = 513.5058
    peptide = Peptide('KAY', 'H', 'Dbz')

    assert abs(true_mass - peptide.mass) <= .01

    # H, Nbz
    true_mass = 539.5878
    peptide = Peptide('KAY', 'H', 'Nbz')

    assert abs(true_mass - peptide.mass) <= .01

def test_calc_sequence_mass_with_non_canoncials():
    # Short peptide
    true_mass = 1022.4653
    peptide = Peptide('ZXOUB', 'H', 'OH')

    assert abs(true_mass - peptide.mass) <= .01

def test_calc_sequence_mass_with_non_canoncials_and_special_termini():
    # Short peptide
    true_mass = 643.6328
    peptide = Peptide('KAYZ', 'H', 'MPAA')

    assert abs(true_mass - peptide.mass) <= .01