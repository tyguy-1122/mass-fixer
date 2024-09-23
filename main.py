###########################################################################
# MassFixer is a tool that is used to determine potential causes of
# differences between expected and observed masses of synthetic peptides
#
# Created by Sam Scherer and Tyler Jones
###########################################################################

import json
import re
from colorama import Fore, Style, init as colorama_init
from source import delta_finder
from source.peptide import Peptide
from source.delta import Delta
from source.delta_set import DeltaSet

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

    # Non-canonical amino acids or protecting groups?
    residues = json.loads(open('data/residues.json').read())
    non_canonicals = [residue for residue in residues if residues[residue]['non_canonical']]
    while True:
        have_non_canonical = input('Does your sequence include any non-canonical amino acids or protecting groups that are not removed during cleavage? (yes/no): ').lower()
        if have_non_canonical == 'yes':
            print('The following are the pre-defined non-canonical residues that you can reference in your sequence input:')
            print('Name\t\tSymbol\t\tMass')
            for non_canonical in non_canonicals:
                print(residues[non_canonical]['name'] + '\t\t' + non_canonical + '\t\t' + str(residues[non_canonical]['mass']))
            non_canonical_num = 1
            while True:
                have_undefined_non_canoncials = input('Do you have any residues that are not included in this list? (yes/no)').lower()
                # Allow user to define new non-canonical residues
                if have_undefined_non_canoncials == 'yes':
                    new_non_canoncical_to_define = True
                    while new_non_canoncical_to_define:
                        new_non_canonical = input('Enter the name of your non-canonical residue')
                        while True:
                            new_non_canoncial_mass = input('Enter the mass of your non-canonical')
                            try:
                                new_non_canoncial_mass = float(new_non_canoncial_mass)
                                break
                            except(Exception):
                                print('The mass for the non-canonical residue must be a number. Please try again.')
                        residues[new_non_canonical] = {'symbol': str(non_canonical_num), 'mass': new_non_canoncial_mass}
                        while True:
                            more_to_define = input('Do you have more non_canonicals to define? (yes/no)').lower()
                            if more_to_define == 'yes':
                                new_non_canoncical_to_define = True
                                non_canonical_num += 1
                            elif more_to_define == 'no':
                                new_non_canoncical_to_define = False
                                print('Here is updated list of non-canonical residues you can reference in your sequence input: ')
                                print('Name\t\tSymbol\t\tMass')
                                for non_canonical in non_canonicals:
                                    print(residues[non_canonical]['name'] + '\t\t' + non_canonical + '\t\t' + str(residues[non_canonical]['mass']))
                                break
                            print('Please answer yes/no. Try again.')
                        break
                    break
                elif have_undefined_non_canoncials == 'no':
                    break
                print('Please answer \"yes\" or \"no\".')
            break
        elif have_non_canonical == 'no':
            break
        print('Please answer \"yes\" or \"no\".')


    # Get sequence
    while True:
        sequence = input('Enter the desired sequence of your peptide: ').upper()
        possible_residue_symbols = ''.join(residues.keys())
        if re.search(f"[^{possible_residue_symbols}]", sequence):
            print('You have entered an invalid sequence. Please try again.')
            continue
        break

    # Get N-Terminus
    termini_species_masses = json.loads(open('data/termini_species.json').read())
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

    return sequence, observed_mass, confidence, n_terminus, c_terminus, non_canonicals

def print_intro():
    '''
    Prints ascii art and a short description of the use cases for the tool
    '''
    print('#############################################################')
    print('#############################################################')
    print(open('resources/mass-fixer-ascii-art.txt').read())
    print('MassFixer is a tool for detecting possible reasons for differences between')
    print('observed and expected masses of synthetic peptides. It is currently capable')
    print('of handling the following potential synthesis issues: \n - residue deletions\n - truncations\n')
    print('#############################################################')
    print('#############################################################\n')


if __name__ == '__main__':
    colorama_init()
    print_intro()

    sequence, obs_mass, uncertainty, n_terminus, c_terminus, non_canonicals = get_user_input()

    # Uncomment this block for debugging and comment out get_user_input
    # obs_mass = 2
    # sequence = 'AAGLT'
    # uncertainty = 2

    # Get masses of residues
    # aa_masses = json.loads(open('data/aa_masses.json').read())
    # termini_species_masses = json.loads(open('data/termini_species_masses.json').read())
    # for non_canonical in non_canonicals:
    #     aa_masses[non_canonicals[non_canonical]['symbol']] = non_canonicals[non_canonical]['mass']
    
    peptide = Peptide(sequence, n_terminus, c_terminus)
    target_mass = peptide.mass - obs_mass
    
    # Validate expected mass
    while True:
        ans = input(f'We calculated your expected mass to be {peptide.mass} and your delta to be {target_mass}. Does that look correct? (yes/no): ').lower()
        if ans == 'no':
            print("Check your math! ...or fix the bug in this program. Goodbye!")
            exit()
        elif ans == 'yes':
            break
        print('Please answer \"yes\" or \"no\".')

    # Find solutions
    print('Calculating...')
    solutions = delta_finder.get_solutions(peptide, target_mass, uncertainty)

    # Show solutions to user
    for solution in solutions:
        print('\n')
        print(solution)
        print('\n')
    # print_pretty_solutions(solutions, sequence, expected_mass, n_terminus, c_terminus)