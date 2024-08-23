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

def calc_expected_mass(sequence: str, aa_masses: dict):
    '''
    calculates expected peptide mass from input sequence
    '''
    total_mass = 0
    for aa in sequence:
        total_mass += aa_masses[aa]
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

def print_pretty_solutions(ugly_solutions, sequence, expected_mass):
    '''
    Prints a solution in a table of sequences (where deleted residues are shown in red), and expected mass
    '''
    if not ugly_solutions:
        print(f'{Fore.RED} No potential solutions found, try increasing your confidence interval.')
        exit()
    print(f'{Fore.BLUE}Sequence{' ' * (len(sequence) - len('Sequence ') if len(sequence) > len('Sequence ') else 0)}\tExpected Mass{Style.RESET_ALL}')
    for ugly_solution in ugly_solutions:
        deleted_residues = ugly_solution['residues']
        pretty_sequence = ''
        for i in range(len(sequence)):
            if i in deleted_residues:
                pretty_sequence += f'{Fore.RED}{sequence[i]}{Style.RESET_ALL}'
            else:
                pretty_sequence += sequence[i]
        print(pretty_sequence + '\t' + str(expected_mass - ugly_solution['mass']))
        

def get_user_input():
    print('#############################################################')
    print('#############################################################')
    print(open('mass-fixer-ascii-art.txt').read())
    print('MassFixer is a tool for detecting possible reasons for differences between')
    print('observed and expected masses of synthetic peptides. It is currently capable')
    print('of handling the following potential synthesis issues: \n - residue deletions\n - truncations\n')
    print('#############################################################')
    print('#############################################################\n')
    
    # Get sequence
    while True:
        sequence = input('Enter the desired sequence of your peptide: ').upper()
        if re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", sequence):
            print('You have entered an invalid sequence. Please try again.')
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

    return sequence, observed_mass, confidence


if __name__ == '__main__':
    # Enable colored terminal output
    colorama_init()

    # Get amino acid masses
    with open ('aa_masses.json', 'r') as aa_json:
        aa_masses = json.load(aa_json)

    # Get user input
    sequence, obs_mass, uncertainty = get_user_input()
    # Uncomment this block for debugging and comment out get_user_input
    # obs_mass = 243
    # sequence = 'AAGLT'
    # uncertainty = 2
    
    # Calc expected mass
    expected_mass = calc_expected_mass(sequence, aa_masses)
    delta_mass = expected_mass - obs_mass

    ans = input(f'We calculated your expected mass to be {expected_mass} and your delta to be {delta_mass}. Does that look correct? (yes/no): ').lower()
    if ans == 'no':
        print("Check your math! ...or fix the bug in this program. Goodbye!")
        exit()

    # Find solutions
    print('Calculating...')
    solutions = delta_mass_combinations(sequence, delta_mass, aa_masses)

    # Show solutions to user
    print_pretty_solutions(solutions, sequence, expected_mass)