def calc_sequence_mass(sequence: str, n_terminus: str, c_terminus: str, aa_masses: dict, termini_species_masses: dict):
    '''
    calculates expected peptide mass from input sequence
    '''
    total_mass = 0
    for aa in sequence:
        total_mass += aa_masses[aa]
    
    total_mass += termini_species_masses[n_terminus]
    total_mass += termini_species_masses[c_terminus]

    return total_mass