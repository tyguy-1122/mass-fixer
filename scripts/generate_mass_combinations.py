import json
from collections import deque

residue_masses = json.loads(open('data/residues.json').read())
# residue_masses = {'A': 5, 'B': 2, 'C': 1}
residues = [residue for residue in residue_masses]
residues_sorted = sorted(residues, key=lambda x: residue_masses[x], reverse=True)

deltas = []
root = {'mass': 0, 'residues': [], 'position': 0}
queue = deque()
queue.append(root)
max_deletions = 4

while queue:
    curr_node = queue.popleft()

    # Condition for solution
    deltas.append(curr_node)
    
    # Termination condtion
    if len(curr_node['residues']) >= max_deletions:
        continue
    
    # Continue recursive search
    for i in range (curr_node['position'], len(residues_sorted)):
        next_mass = residue_masses[residues_sorted[i]]
        queue.append({'mass': curr_node['mass'] + next_mass, 'residues': curr_node['residues'] + [residues_sorted[i]], 'position': i})

print(len(deltas))
