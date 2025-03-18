import pyrosetta
from Bio.PDB import PDBIO, Structure
from pyrosetta import pose_from_pdb, get_fa_scorefxn

def compute_rosetta_energy(structure: Structure.Structure, output_pdb="temp.pdb"):
    """
    Compute the Rosetta energy score for a given Biopython structure.
    
    Args:
        structure (Bio.PDB.Structure.Structure): The input protein structure.
        output_pdb (str): Temporary PDB filename for conversion to PyRosetta.

    Returns:
        float: The computed Rosetta energy score.
    """
    # Initialize PyRosetta
    pyrosetta.init()

    # Save Biopython structure to a temporary PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

    # Load the structure into PyRosetta
    pose = pose_from_pdb(output_pdb)

    # Get the full-atom Rosetta scoring function
    scorefxn = get_fa_scorefxn()

    # Compute and return the energy score
    energy = scorefxn(pose)
    return energy

