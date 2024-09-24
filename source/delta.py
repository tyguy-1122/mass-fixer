from enum import Enum

class DeltaType(Enum):
    DELETION = "DELETION"
    TRUNCATION = "TRUNCATION"

class Delta:
    def __init__(self, mass, type: DeltaType, likelihood: int, description=""):
        """
        Initialize a Delta instance with mass, type, and description.
        """
        if not isinstance(type, DeltaType):
            raise ValueError("Type must be an instance of DeltaType Enum")
        
        self.mass = mass
        self.type = type
        self.description = description
        self.likelihood = likelihood

    def __repr__(self):
        return f"Delta(mass={self.mass}, type={self.type.value}, description='{self.description}')"

    def __hash__(self):
        return hash(self.description)
    
    def __eq__(self, other):
        return self.mass == other.mass and self.type == other.type and self.description == other.description