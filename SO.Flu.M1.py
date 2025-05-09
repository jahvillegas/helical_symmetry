import sys
import os
import random
import math
import csv
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue
import numpy as np
from core.helical_symmetry import compute_helical_parameters  # Make sure this import is correct
from core.energy import *  # Make sure this import is correct

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
    #pdb_file = "src/data/R3K_16_AU_centered.pdb"  # Replace with your input PDB file path
    pdb_file = "data/7jm3.pdb"  # Replace with your input PDB file path


    pitch = 35
    n_units = 21

    rise = pitch / n_units  
    twist = 360 / n_units
    num_units = 20  # Number of asymmetric units
    best_energy = 10000000000000000000000
    displacement = 0

    dof = [rise, twist, displacement]
    dr_max = [1, 5, 1]
    #dr_max = [0, 0, 0]

    #f = open('data/monte_carlo.csv', 'w')
    
    # Define the ranges for the parameters
    pitch_range = (20, 45)      
    num_units_range = (10, 25)  
    displacement_range = (-10, 10)

    n_steps = 10       
    
    for i in range(n_steps):
        # Randomly select parameters within the defined ranges
        pitch = random.uniform(pitch_range[0], pitch_range[1])
        num_units = random.randint(num_units_range[0], num_units_range[1])
        displacement = random.uniform(displacement_range[0], displacement_range[1])

        # Load asymmetric unit
        structure = load_pdb(pdb_file)

        #Center protein with rotation asix at origin
        translation_vector = np.array([-200, -200, -200])
        structure = displace_protein(structure, translation_vector)

        #Displace protein from center of mass to origin by the given distance
        structure = move_toward_origin(structure, 0)
        # Create helical assembly
        assembly = create_helical_assembly(structure, twist, rise, num_units)

        # Save the helical assembly
        save_pdb(assembly, "data/M1_assembly_2.pdb")

        # Compute Rosetta energy
        energy_m = compute_rosetta_energy(assembly, minimize=False)
        print(f"Iteration {i}: Energy_m = {energy_m:.2f}, pitch = {pitch:.2f}, num_units = {num_units:.2f}, displacement = {displacement:.2f}")
        if energy_m < best_energy:
            best_energy = energy_m
            best_assembly = assembly
            print(f"New best energy: {best_energy:.2f}")
        selected_dof = random.randint(0, len(dof) - 1)
        dof_old = dof[selected_dof]
        change =  random.uniform(-1 * dr_max[selected_dof], dr_max[selected_dof])
        dof[selected_dof] = dof[selected_dof] + change
    
        # Load asymmetric unit
        structure = load_pdb(pdb_file)
    
        #Center protein with rotation asix at origin
        translation_vector = np.array([-200, -200, -200])
        structure = displace_protein(structure, translation_vector)
    
        #Displace protein from center of mass to origin by the given distance
        structure = move_toward_origin(structure, dof[2])
        # Create helical assembly
        assembly = create_helical_assembly(structure, dof[1], dof[0], num_units)
    
        # Save the helical assembly
        save_pdb(assembly, "data/M1_assembly_4.pdb")

        energy_n = compute_rosetta_energy(assembly, minimize=False)
        delta = energy_n - energy_m
        if delta > 0:
            zeta = random.uniform(0,1)
            metropolis = math.exp((-1 * delta) > zeta)
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

        #writer = csv.writer(f)
        #writer.writerow([i, energy_m, dof])
    
    print(f"Rosetta Energy Score: {energy_m:.2f}")
if __name__ == "__main__":
    main()
 

        
        
 
        