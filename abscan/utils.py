import os
import sys
import subprocess
import logging
import numpy as np

logger = logging.getLogger("abscan.utils")

def get_obabel_path():
    """Dynamically locate the obabel executable in the current environment."""
    env_dir = os.path.dirname(sys.executable)
    obabel_path = os.path.join(env_dir, "obabel")
    if os.path.exists(obabel_path):
        return obabel_path
    
    # Fallback to system path
    import shutil
    system_obabel = shutil.which("obabel")
    if system_obabel:
        return system_obabel
        
    raise RuntimeError("obabel executable not found in conda environment or system PATH.")

def get_ligand_coords(ligand_pdb_path):
    """
    Parse a ligand PDB file and return a numpy array of coordinates of all atoms.
    """
    coords = []
    with open(ligand_pdb_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append([x, y, z])
                except ValueError:
                    continue
    if not coords:
        raise ValueError(f"No atoms found in ligand PDB file: {ligand_pdb_path}")
    return np.array(coords)

def get_vina_box(coords, padding=10.0):
    """
    Calculate the center and size of the box enclosing the ligand.
    """
    center = np.mean(coords, axis=0).tolist()
    box_size = (np.max(coords, axis=0) - np.min(coords, axis=0) + padding).tolist()
    return center, box_size

def run_command(cmd_list, description=""):
    """Run a system command and log errors if it fails."""
    logger.debug(f"Running command: {' '.join(cmd_list)}")
    result = subprocess.run(cmd_list, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Error during {description or 'command'}: {result.stderr}")
        raise RuntimeError(f"Command failed: {' '.join(cmd_list)}\nError: {result.stderr}")
    return result.stdout

def prepare_receptor(input_pdb, output_pdbqt):
    """
    Convert receptor PDB to PDBQT format using obabel.
    Uses '-xr' to treat as rigid receptor and adds gasteiger charges.
    """
    obabel = get_obabel_path()
    cmd = [
        obabel,
        input_pdb,
        "-O", output_pdbqt,
        "-xr",
        "-p", "7.4",
        "--partialcharge", "gasteiger"
    ]
    run_command(cmd, f"receptor preparation ({input_pdb} -> {output_pdbqt})")

def prepare_ligand(input_pdb, output_pdbqt):
    """
    Convert ligand PDB to PDBQT format using obabel.
    Detects torsions automatically and adds gasteiger charges.
    """
    obabel = get_obabel_path()
    cmd = [
        obabel,
        input_pdb,
        "-O", output_pdbqt,
        "-p", "7.4",
        "--partialcharge", "gasteiger"
    ]
    run_command(cmd, f"ligand preparation ({input_pdb} -> {output_pdbqt})")
