###########################################################################
# MassFixer is a tool that is used to determine potential causes of
# differences between expected and observed masses of synthetic peptides
#
# Created by Sam Scherer and Tyler Jones
###########################################################################

import json
import re
from collections import deque
from colorama import Fore, Style, init as colorama_init

def calc_expected_mass(sequence: str, n_terminus: str, c_terminus: str, aa_masses: dict, termini_species_masses: dict):
    '''
    calculates expected peptide mass from input sequence
    '''
    total_mass = 0
    for aa in sequence:
        total_mass += aa_masses[aa]
    
    total_mass += termini_species_masses[n_terminus]
    total_mass += termini_species_masses[c_terminus]

    return total_mass

def delta_mass_combinations(sequence: str, delta_mass: float, aa_masses: dict):
    '''
    Returns list of possible combinations of deletions that would result in the delta mass.
    
    Uses a breadth-first-search algorithm to search the solution space.
    '''
    solutions = []
    root = {'mass': 0, 'residues': [], 'position': -1}
    queue = deque()
    queue.append(root)

    # Perform BFS
    while queue:
        curr_node = queue.popleft()

        # Condition for solution
        if delta_mass - uncertainty < curr_node['mass'] and curr_node['mass'] < delta_mass + uncertainty:
            solutions.append(curr_node)
        
        # Termination condtion
        elif curr_node['mass'] >= uncertainty + delta_mass:
            continue
        
        # Continue recursive search
        for i in range (curr_node['position'] + 1, len(sequence)):
            next_mass = aa_masses[sequence[i]]
            queue.append({'mass': curr_node['mass'] + next_mass, 'residues': curr_node['residues'] + [i], 'position': i})
    
    return solutions

def print_pretty_solutions(ugly_solutions, sequence, expected_mass, n_terminus, c_terminus):
    '''
    Prints a solution in a table of sequences and expected masses (where deleted residues are shown in red).
    '''
    if not ugly_solutions:
        print(f'{Fore.RED} No potential solutions found, try increasing your confidence interval.')
        exit()

    sequence_column_length = len(sequence) + len(n_terminus) + len(c_terminus)
    print(f'{Fore.BLUE}Sequence{' ' * (sequence_column_length - len('Sequence ') if sequence_column_length > len('Sequence ') else 0)}\tExpected Mass{Style.RESET_ALL}')
    
    # Loop through the sequence predictions and color deleted residues red
    for ugly_solution in ugly_solutions:
        deleted_residues = ugly_solution['residues']
        pretty_sequence = n_terminus + '--'
        for i in range(len(sequence)):
            if i in deleted_residues:
                pretty_sequence += f'{Fore.RED}{sequence[i]}{Style.RESET_ALL}'
            else:
                pretty_sequence += sequence[i]
        pretty_sequence += '--'
        pretty_sequence += c_terminus
        print(pretty_sequence + '\t' + str(expected_mass - ugly_solution['mass']))
        

def get_user_input():
    '''
    Gets user input about the peptide and observed mass
    '''
    print('#############################################################')
    print('#############################################################')
    print(open('mass-fixer-ascii-art.txt').read())
    print('MassFixer is a tool for detecting possible reasons for differences between')
    print('observed and expected masses of synthetic peptides. It is currently capable')
    print('of handling the following potential synthesis issues: \n - residue deletions\n - truncations\n')
    print('#############################################################')
    print('#############################################################\n')

    # Non-canonical amino acids or protecting groups?
    while True:
        non_canonical = input('Does your sequence include any non-canonical amino acids or protecting groups that are not removed during cleavage? (yes/no): ').lower()
        if ans == 'yes':
            print('TODO: Implement non-canonicals')
            break
        elif ans == 'no':
            break
        print('Please answer \"yes\" or \"no\".')


    # Get sequence
    while True:
        sequence = input('Enter the desired sequence of your peptide: ').upper()
        if re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", sequence):
            print('You have entered an invalid sequence. Please try again.')
            continue
        break

    # Get N-Terminus
    termini_species_masses = json.loads(open('termini_species_masses.json').read())
    print('Here are the possible chemical species for your termini for reference: ')
    print(list(termini_species_masses.keys()))
    while True:
        n_terminus = input('Enter the chemical species at the N-terminus of your peptide: ')
        if not n_terminus in termini_species_masses:
            print('You have entered an invalid n_terminus. Please try again.')
            continue
        break

    # Get C-Terminus
    while True:
        c_terminus = input('Enter the chemical species at the C-terminus of your peptide: ')
        if not c_terminus in termini_species_masses:
            print('You have entered an invalid c_terminus. Please try again.')
            continue
        break

    # Get observed mass
    while True:
        observed_mass = input('Enter the observed mass of your peptide: ')
        try:
            observed_mass = float(observed_mass)
            break
        except(Exception):
            print('Observed mass must be a number. Please try again.')

    # Get uncertainty
    while True:
        confidence = input('Enter the confidence of your observed mass in Daltons: ')
        try:
            confidence = float(confidence)
            break
        except(Exception):
            print('Confidence must be a number. Please try again.')

    return sequence, observed_mass, confidence, n_terminus, c_terminus


if __name__ == '__main__':
    # Enable colored terminal output
    colorama_init()

    # Get masses
    aa_masses = json.loads(open('aa_masses.json').read())
    termini_species_masses = json.loads(open('termini_species_masses.json').read())

    # Get user input
    sequence, obs_mass, uncertainty, n_terminus, c_terminus = get_user_input()
    # Uncomment this block for debugging and comment out get_user_input
    # obs_mass = 2
    # sequence = 'AAGLT'
    # uncertainty = 2
    
    # Calc expected mass
    expected_mass = calc_expected_mass(sequence, n_terminus, c_terminus, aa_masses, termini_species_masses)
    delta_mass = expected_mass - obs_mass

    while True:
        ans = input(f'We calculated your expected mass to be {expected_mass} and your delta to be {delta_mass}. Does that look correct? (yes/no): ').lower()
        if ans == 'no':
            print("Check your math! ...or fix the bug in this program. Goodbye!")
            exit()
        elif ans == 'yes':
            break
        print('Please answer \"yes\" or \"no\".')

    # Find solutions
    print('Calculating...')
    solutions = delta_mass_combinations(sequence, delta_mass, aa_masses)

    # Show solutions to user
    print_pretty_solutions(solutions, sequence, expected_mass, n_terminus, c_terminus)