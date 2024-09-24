from source.delta import Delta

class DeltaSet:
    def __init__(self, deltas: frozenset[Delta], mass: float = None, likelihood: int = None):
        self.deltas = deltas
        self.mass = mass if mass else self.calculate_mass()
        self.likelihood = likelihood if likelihood else self.calculate_likelihood()
    
    def calculate_mass(self):
        total_mass = 0
        for delta in self.deltas:
            total_mass += delta.mass

        return total_mass
    
    def calculate_likelihood(self):
        likelihood = 0
        for delta in self.deltas:
            likelihood += delta.likelihood
        
        return likelihood

    def __hash__(self):
        return hash(self.deltas)
    
    def __eq__(self, other):
        return hash(self.deltas) == hash(other.deltas) and self.mass == other.mass  
    
    def __repr__(self):
        return f"DeltaSet(mass={self.mass}, deltas='{self.deltas}')"