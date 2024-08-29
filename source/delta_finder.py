from collections import deque

def find_possible_deltas(sequence: str, delta_mass: float, residue_masses: dict, confidence: float):
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
        if delta_mass - confidence < curr_node['mass'] and curr_node['mass'] < delta_mass + confidence:
            solutions.append(curr_node)
        
        # Termination condtion
        elif curr_node['mass'] >= confidence + delta_mass:
            continue
        
        # Continue recursive search
        for i in range (curr_node['position'] + 1, len(sequence)):
            next_mass = residue_masses[sequence[i]]
            queue.append({'mass': curr_node['mass'] + next_mass, 'residues': curr_node['residues'] + [i], 'position': i})
    
    return solutions