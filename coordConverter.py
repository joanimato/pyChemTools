from typing import List, Tuple
from dataclasses import dataclass
from numpy.typing import ArrayLike
import numpy as np

@dataclass
class Coordinates:
    cartesianCoords: ArrayLike 
    internalCoords: List[str]
    __COVALENT_RADII = {
    'H': 0.31,  # Hydrogen
    'C': 0.76,  # Carbon
    'O': 0.66,  # Oxygen
    'N': 0.71,  # Nitrogen
}

    @staticmethod
    def create_adjacency_list(cartesianCoords, atom_types, tolerance=0.45):
        """
        Create an adjacency matrix from Cartesian coordinates and atomic types.
        
        Parameters:
        - cartesian_coords: List of lists of Cartesian coordinates [[x1, y1, z1], [x2, y2, z2], ...].
        - atom_types: List of atomic symbols corresponding to the Cartesian coordinates ['C', 'O', 'H', ...].
        - tolerance: A small margin added to the sum of covalent radii to allow flexibility in bonding.
        
        Returns:
        - adjacency_list
        """
        num_atoms = len(cartesianCoords)
        adjacency_list = []

        # Loop over all pairs of atoms to calculate distances and check for bonds
        for i in range(num_atoms):
            for j in range(i+1, num_atoms):
                # Get covalent radii for the atom types
                radius_i = Coordinates.__COVALENT_RADII.get(atom_types[i], 0)
                radius_j = Coordinates.__COVALENT_RADII.get(atom_types[j], 0)

                # Calculate distance between atom i and atom j
                dist = Coordinates.__distance(cartesianCoords[i], cartesianCoords[j])

                # If the distance is less than or equal to the sum of covalent radii + tolerance, they are bonded
                if dist <= (radius_i + radius_j + tolerance):
                    adjacency_list.append((i, j))

        return adjacency_list
    

    @staticmethod
    def __distance(x,y):
        return np.linalg.norm(np.array(x) - np.array(y))
    
    @staticmethod
    def __angle(atom1, atom2, atom3):
        """Calculate bond angle between three atoms."""
        vec1 = np.array(atom1) - np.array(atom2)
        vec2 = np.array(atom3) - np.array(atom2)
        
        cosine_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        angle = np.arccos(cosine_angle)
        
        return np.degrees(angle)
    
    @staticmethod
    def __dihedral(atom1, atom2, atom3, atom4):
        """Calculate dihedral (torsion) angle between four atoms."""
        vec1 = np.array(atom2) - np.array(atom1)
        vec2 = np.array(atom3) - np.array(atom2)
        vec3 = np.array(atom4) - np.array(atom3)
        
        norm1 = np.cross(vec1, vec2)
        norm2 = np.cross(vec2, vec3)
        
        norm1 /= np.linalg.norm(norm1)
        norm2 /= np.linalg.norm(norm2)
        
        m1 = np.cross(norm1, vec2 / np.linalg.norm(vec2))
        
        x = np.dot(norm1, norm2)
        y = np.dot(m1, norm2)
        
        return np.degrees(np.arctan2(y, x))

    @staticmethod
    def generate_internal_coordinates(cartesian_coords, adjacency_list):
        """
        Generate internal coordinates for input into computational chemistry software.
        Format: 
        Atom Index 
        Connected Atom 1 bond_length 1
        Connected Atom 2 bond_length 1 angle 2
        Connected Atom 3 bond_length 1 angle 2 dihedral 3
        """
        internal_coords = []
        
        for i in range(len(cartesian_coords)):
            if i == 0:
                # First atom, no internal coordinates (typically origin)
                internal_coords.append(f"{i + 1}")
            elif i == 1:
                # Second atom, has one bond
                j = adjacency_list[0][0] if adjacency_list[0][1] == i else adjacency_list[0][1]
                bond_length = Coordinates.__distance(cartesian_coords[i], cartesian_coords[j])
                internal_coords.append(f"{i + 1} {j + 1} {bond_length:.6f}")
            elif i == 2:
                # Third atom, has one bond and one angle
                bond_pair = adjacency_list[i - 1]
                j = bond_pair[0] if bond_pair[1] == i else bond_pair[1]
                bond_length = Coordinates.__distance(cartesian_coords[i], cartesian_coords[j])
                
                # Angle calculation
                k = adjacency_list[0][0] if adjacency_list[0][1] == j else adjacency_list[0][1]
                bond_angle = Coordinates.__angle(cartesian_coords[k], cartesian_coords[j], cartesian_coords[i])
                internal_coords.append(f"{i + 1} {j + 1} {bond_length:.6f} {k + 1} {bond_angle:.6f}")
            else:
                # Fourth and beyond, bond, angle, and dihedral
                bond_pair = adjacency_list[i - 1]
                j = bond_pair[0] if bond_pair[1] == i else bond_pair[1]
                bond_length = Coordinates.__distance(cartesian_coords[i], cartesian_coords[j])
                
                # Angle calculation
                k = adjacency_list[j - 1][0] if adjacency_list[j - 1][1] == j else adjacency_list[j - 1][1]
                bond_angle = Coordinates.__angle(cartesian_coords[k], cartesian_coords[j], cartesian_coords[i])
                
                # Dihedral calculation
                l = adjacency_list[k - 1][0] if adjacency_list[k - 1][1] == k else adjacency_list[k - 1][1]
                torsion_angle = Coordinates.__dihedral(cartesian_coords[l], cartesian_coords[k], cartesian_coords[j], cartesian_coords[i])
                
                internal_coords.append(f"{i + 1} {j + 1} {bond_length:.6f} {k + 1} {bond_angle:.6f} {l + 1} {torsion_angle:.6f}")
        
        return internal_coords

    @staticmethod
    def fromCartesian(coords: ArrayLike, atom_types):
        """Generate a Coordinates object from an array of cartesian coordinates"""
        cartesian = coords 
        adjMatrix = Coordinates.create_adjacency_list(coords, atom_types)
        internals = Coordinates.generate_internal_coordinates(cartesian, adjMatrix)
        return Coordinates(cartesian, internals)



cartesian_coords = [
    [0.0, 0.0, 0.0],   # Atom 1 (C1)
    [1.54, 0.0, 0.0],  # Atom 2 (C2)
    [0.0, 1.0, 0.0],   # Atom 3 (H1)
    [0.0, -1.0, 0.0],  # Atom 4 (H2)
    [1.0, 0.0, 1.0],   # Atom 5 (H3)
    [2.54, 1.0, 0.0],  # Atom 6 (H4)
    [2.54, -1.0, 0.0], # Atom 7 (H5)
    [1.54, 0.0, -1.0]  # Atom 8 (H6)
]

atom_types = ['C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']

adjacency_list = Coordinates.create_adjacency_list(cartesian_coords, atom_types)
internal_coords = Coordinates.generate_internal_coordinates(cartesian_coords, atom_types, adjacency_list)

coords = Coordinates.fromCartesian(cartesian_coords, atom_types)

for i in internal_coords:
    print(i)