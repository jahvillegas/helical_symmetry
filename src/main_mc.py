import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue
import numpy as np
from core.helical_symmetry import compute_helical_parameters  # Make sure this import is correct
from core.energy import *  # Make sure this import is correct

import random
import math
import csv


# 🟢 Load PDB File
def load_pdb(file_path):
    parser = PDBParser()
    structure = parser.get_structure("asymmetric_unit", file_path)
    return structure

# 🟢 Create Helical Assembly
def create_helical_assembly(structure, twist, rise, num_units):
    rotation_matrix, translation_vector = compute_helical_parameters(twist, rise)
    model = structure[0]
    assembly = Structure.Structure("helical_assembly")
    new_model = Model.Model(0)
    atom_serial = 1  # Track unique atom serial numbers

    for i in range(num_units):
        new_chain = Chain.Chain(chr(65 + i))  # Use A, B, C,... for chain IDs
        for residue in model.get_residues():
            # 🟢 Create a unique residue for each copy
            new_residue = Residue.Residue(
                (" ", residue.get_id()[1] + i * 100, " "),  # Unique residue ID
                residue.get_resname(),
                residue.get_segid()
            )
            for atom in residue.get_atoms():
                # Transform coordinates
                coord = atom.get_coord()
                transformed_coord = np.dot(np.linalg.matrix_power(rotation_matrix, i), coord) + i * translation_vector

                # 🟢 Create a new atom with a unique serial number
                new_atom = atom.copy()
                new_atom.set_coord(transformed_coord)
                new_atom.serial_number = atom_serial  # Unique serial number
                atom_serial += 1  # Increment serial number

                # Add transformed atom to new residue
                new_residue.add(new_atom)

            # Add new residue to chain
            new_chain.add(new_residue)

        # Add the new chain to the model
        new_model.add(new_chain)

    # Add the model to the assembly
    assembly.add(new_model)
    return assembly

# 🟢 Save Transformed PDB
def save_pdb(structure, output_file):
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)
    print(f"Helical assembly saved as {output_file}")

def displace_protein(structure, displacement_vector):
    """
    Displace all atoms in a protein structure by a given translation vector.
    
    Args:
        structure (Bio.PDB.Structure.Structure): The protein structure to translate.
        displacement_vector (array-like): A 3D translation vector [dx, dy, dz] in Ångströms.

    Returns:
        Bio.PDB.Structure.Structure: Translated protein structure.
    """
    for atom in structure.get_atoms():
        new_coord = atom.get_coord() + displacement_vector
        atom.set_coord(new_coord)
    
    return structure

def move_toward_origin(structure, distance: float):
    """
    Move a protein structure toward the origin by a specified distance.

    Args:
        structure (Bio.PDB.Structure.Structure): The protein structure to move.
        distance (float): The distance to move toward the origin (in Ångströms).

    Returns:
        Bio.PDB.Structure.Structure: Translated protein structure.
    """
    # 🟢 Compute center of mass (COM)
    atom_coords = np.array([atom.get_coord() for atom in structure.get_atoms()])
    center_of_mass = atom_coords.mean(axis=0)

    # 🟢 Compute unit direction vector from COM to the origin
    direction = -center_of_mass / np.linalg.norm(center_of_mass)

    # 🟢 Compute translation vector
    translation_vector = direction * distance

    # 🟢 Apply translation
    for atom in structure.get_atoms():
        new_coord = atom.get_coord() + translation_vector
        atom.set_coord(new_coord)

    return structure

# 🟢 Main Function
def main():
    #pdb_file = "src/data/M1_fibril.pdb"  # Replace with your input PDB file path
    pdb_file = "src/data/M1_fibril.pdb"  # Replace with your input PDB file path

    #rise = 4.155  
    #twist = 22.830    
    #displacement = 0

    rise = 5.972 
    twist = 32.818    
    displacement = 15

    n_steps = 20
    best_energy = 10000000000000000000000

    dof = [rise, twist, displacement]
    dr_max = [1, 5, 1]
    #dr_max = [0, 0, 0]

    with open('src/data/monte_carlo.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Step", "Energy", "DOF"])

    for i in range(n_steps):
        ...
        writer.writerow([i, energy_m, dof])


    num_units = 9  # Number of asymmetric units

    for i in range(n_steps): 

    
        # Load asymmetric unit
        structure = load_pdb(pdb_file)
    
        #Center protein with rotation asix at origin
        translation_vector = np.array([-336.4, -336.4, -336.4])
        structure = displace_protein(structure, translation_vector)
    
        #Displace protein from center of mass to origin by the given distance
        structure = move_toward_origin(structure, dof[2])
        # Create helical assembly
        assembly = create_helical_assembly(structure, dof[1], dof[0], num_units)
    
        # Save the helical assembly
        save_pdb(assembly, "src/data/R3K_N16_helical_assembly_4.pdb")
    
        # Compute Rosetta energy
        energy_m  = compute_rosetta_energy(assembly)

        selected_dof = random.randint(0, len(dof) - 1)
        dof_old = dof[selected_dof]
        change =  random.uniform(-1 * dr_max[selected_dof], dr_max[selected_dof])
        dof[selected_dof] = dof[selected_dof] + change
    
        # Load asymmetric unit
        structure = load_pdb(pdb_file)
    
        #Center protein with rotation asix at origin
        translation_vector = np.array([-336.4, -336.4, -336.4])
        structure = displace_protein(structure, translation_vector)
    
        #Displace protein from center of mass to origin by the given distance
        structure = move_toward_origin(structure, dof[2])
        # Create helical assembly
        assembly = create_helical_assembly(structure, dof[1], dof[0], num_units)
    
        # Save the helical assembly
        save_pdb(assembly, "src/data/R3K_N16_helical_assembly_4.pdb")

        energy_n = compute_rosetta_energy(assembly)
        delta = energy_n - energy_m
        
 
        if delta > 0:
            zeta = random.uniform(0,1)
            metropolis = math.exp((-1 * delta) / Temp ) > zeta
            print(selected_dof)
            #dof[selected_dof] = dof_old
            if(metropolis == False):
               print(f"Rosetta Energy Score: {energy_m:.2f}")
               print(f"Rosetta Energy Score: {energy_n:.2f}")
               print(dof_old)
               print(dof)
               dof[selected_dof] = dof_old

        print(f"Rosetta Energy Score: {energy_m:.2f}")
        print(dof)

        writer = csv.writer(f)
        writer.writerow([i, energy_m, dof])

    f.close()
if __name__ == "__main__":
    main()
