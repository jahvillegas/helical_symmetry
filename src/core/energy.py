import pyrosetta
from Bio.PDB import PDBIO, Structure
from pyrosetta import pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.optimization import AtomTreeMinimizer
from pyrosetta.rosetta.core.scoring import ScoreFunction
from pyrosetta.rosetta.core.kinematics import MoveMap

def compute_rosetta_energy(structure: Structure.Structure, minimize=False, output_pdb="temp.pdb"):
    """
    Compute the Rosetta energy score for a given Biopython structure, with optional minimization (no repacking).

    Args:
        structure (Bio.PDB.Structure.Structure): The input protein structure.
        minimize (bool): Whether to perform energy minimization before scoring (without repacking).
        output_pdb (str): Temporary PDB filename for conversion to PyRosetta.

    Returns:
        float: The Rosetta energy score (minimized if `minimize=True`).
    """
    # 🟢 Initialize PyRosetta
    pyrosetta.init()

    # 🟢 Save Biopython structure to a temporary PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

    # 🟢 Load the structure into PyRosetta
    pose = pose_from_pdb(output_pdb)

    # 🟢 Get the full-atom Rosetta scoring function
    scorefxn = get_fa_scorefxn()

    # 🟢 Perform energy minimization (NO REPACKING) if requested
    if minimize:
        print("Performing energy minimization (without repacking)...")
        
        # Define a MoveMap (allow backbone minimization but not repacking)
        movemap = MoveMap()
        movemap.set_bb(True)   # Allow backbone minimization
        movemap.set_chi(False) # Disable side-chain movement (NO REPACKING)

        # Setup the minimizer
        minimizer = MinMover()
        minimizer.movemap(movemap)
        minimizer.score_function(scorefxn)
        minimizer.apply(pose)

    # 🟢 Compute and return the energy
    return scorefxn(pose)

