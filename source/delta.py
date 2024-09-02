from enum import Enum

class DeltaType(Enum):
    DELETION = "Deletion"
    TRUNCATION = "Truncation"

class Delta:
    def __init__(self, mass, type: DeltaType, description=""):
        """
        Initialize a Delta instance with mass, type, and description.
        """
        if not isinstance(type, DeltaType):
            raise ValueError("Type must be an instance of DeltaType Enum")
        
        self.mass = mass
        self.type = type
        self.description = description

    def __repr__(self):
        return f"Delta(mass={self.mass}, type={self.type.value}, description='{self.description}')"

    def __hash__(self):
        return hash(self.description)