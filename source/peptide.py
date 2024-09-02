import json

class Peptide:
    termini_species_masses = json.loads(open('data/termini_species_masses.json').read())
    aa_masses = json.loads(open('data/aa_masses.json').read())
    non_canonicals = json.loads(open('data/non-canonical_aas.json').read())
    for non_canonical in non_canonicals:
        aa_masses[non_canonicals[non_canonical]['symbol']] = non_canonicals[non_canonical]['mass']

    def __init__(self, sequence, n_termini_species, c_termini_species, mass=None):
        self.sequence = sequence
        self.n_termini_species = n_termini_species
        self.c_termini_species = c_termini_species
        self.mass = mass if mass else self.calculate_mass()

    def calculate_mass(self):
            '''
            calculates expected peptide mass from input sequence
            '''
            total_mass = 0
            for aa in self.sequence:
                total_mass += self.aa_masses[aa]
            
            total_mass += self.termini_species_masses[self.n_termini_species]
            total_mass += self.termini_species_masses[self.c_termini_species]

            return total_mass