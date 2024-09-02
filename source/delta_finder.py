from collections import deque
import json
from source.peptide import Peptide
from source.delta import Delta, DeltaType
from source.delta_set import DeltaSet

def get_deltas(peptide: Peptide):
    '''
    Gets all of the possible deltas based on the provided sequence
    '''
    # Deletion deltas
    residue_masses = json.loads(open('data/aa_masses.json').read())
    residue_deltas = [Delta(residue_masses[residue], DeltaType.DELETION, str(residue)) for residue in peptide.sequence]

    # TODO: Incomplete Fmoc removal, sodium adducts, incomplete W deprotection, etc.

    return residue_deltas

def generate_delta_sets(deltas: list[Delta], target_mass: float, confidence: float):
    '''
    Generates all possible combinations of deltas for the problem sequence that match the target mass.
    '''
    deltas = sorted(deltas, key=lambda x: x.mass, reverse=True)

    delta_combinations = []
    root = tuple([DeltaSet(frozenset(), 0), 0])
    queue = deque()
    queue.append(root)
    max_deletions = 4
    visited = set()

    while queue:
        curr_node = queue.popleft()

        # Condition for solution
        if (target_mass - confidence <= curr_node[0].mass <= target_mass + confidence):
            delta_combinations.append(curr_node[0])
        
        # Termination condtions
        if len(curr_node[0].deltas) >= max_deletions or curr_node[0].mass >= target_mass + confidence:
            continue
        
        # Continue recursive search
        for i in range (curr_node[1] + 1, len(deltas)):
            deltas_set = curr_node[0].deltas | frozenset([deltas[i]])
            if set(deltas_set) in visited: # Don't do duplicate work
                continue
            visited.add(deltas_set)
            next_mass = deltas[i].mass
            queue.append(tuple([DeltaSet(deltas_set, curr_node[0].mass + next_mass), i]))
    
    return delta_combinations

def get_truncations(peptide: Peptide, target_mass: float, confidence: float):
    '''
    Get all of the N-terminal trucations for the provided peptide without exceeding the target mass.
    '''
    aa_masses = json.loads(open('data/aa_masses.json').read())
    non_canonicals = json.loads(open('data/non-canonical_aas.json').read())
    for non_canonical in non_canonicals:
        aa_masses[non_canonicals[non_canonical]['symbol']] = non_canonicals[non_canonical]['mass']
    truncations = [peptide]
    truncated_mass = 0
    for i in range(1, len(peptide.sequence)):
        truncated_mass += aa_masses[peptide.sequence[i - 1]]
        if truncated_mass > target_mass + confidence:
            break
        truncations.append(Peptide(peptide.sequence[i:len(peptide.sequence)], peptide.n_termini_species, peptide.c_termini_species, truncated_mass))
    return truncations

def get_solutions(peptide: Peptide, target_mass: float, confidence):
    solutions = []
    # TODO: At some point, the most obvious optmization would be to not do the duplicate work of generating the same
    # delta combinations repeatedly for each of the trucations.
    truncations = get_truncations(peptide, target_mass, confidence)

    for truncation in truncations:
        deltas = get_deltas(peptide)
        for delta_combination in generate_delta_sets(deltas, target_mass - (peptide.mass - truncation.mass), confidence):
            solutions.append(DeltaSet(delta_combination.deltas | frozenset([Delta(peptide.mass - truncation.mass, DeltaType.TRUNCATION, f"truncation of \'{peptide.sequence[:len(peptide.sequence) - len(truncation.sequence)]}\'")])))
    
    return solutions