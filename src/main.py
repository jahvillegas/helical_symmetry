from Bio.PDB import PDBParser, Superimposer
import numpy as np
from scipy.optimize import minimize
from src.core.helical_symmetry import compute_helical_parameters


# 🟢 Step 1: Load and Parse PDB File
def load_pdb(file_path):
    parser = PDBParser()
    structure = parser.get_structure("helical_structure", file_path)
    atoms = [atom for atom in structure.get_atoms() if atom.element != 'H']
    return np.array([atom.coord for atom in atoms])


# 🟢 Step 2: Generate Helical Models Based on Parameters
def generate_helical_model(atoms, n1, n2, twist, rise):
    matrix, vector = compute_helical_parameters(n1, n2, twist, rise)
    transformed_atoms = []

    for atom in atoms:
        # Apply helical transformation
        transformed_atom = np.dot(matrix, atom) + vector
        transformed_atoms.append(transformed_atom)

    return np.array(transformed_atoms)


# 🟢 Step 3: Calculate RMSD Between Structures
def calculate_rmsd(coords1, coords2):
    # Align and calculate RMSD using Biopython's Superimposer
    sup = Superimposer()
    sup.set_atoms(coords1, coords2)
    return sup.rms


# 🟢 Step 4: Optimize Helical Parameters to Minimize RMSD
def optimize_helical_parameters(atoms, initial_params):
    def objective(params):
        n1, n2, twist, rise = params
        generated_model = generate_helical_model(atoms, int(n1), int(n2), twist, rise)
        return calculate_rmsd(atoms, generated_model)

    result = minimize(objective, initial_params, method='Nelder-Mead')
    return result.x, result.fun


# 🟢 Main Function
def main():
    pdb_file = "src/data/helical_model.pdb"  # Replace with your PDB file path
    atoms = load_pdb(pdb_file)

    # Initial guesses for n1, n2, twist, rise
    initial_params = [10, 5, 30.0, 4.8]

    # Optimize parameters
    optimized_params, rmsd = optimize_helical_parameters(atoms, initial_params)

    print(f"Optimized Parameters: n1={optimized_params[0]}, n2={optimized_params[1]}, "
          f"twist={optimized_params[2]:.2f}, rise={optimized_params[3]:.2f}")
    print(f"Final RMSD: {rmsd:.2f} Å")


if __name__ == "__main__":
    main()

