# abscan: Alanine Binding Site Mutagenesis Scanning

`abscan` is a modern, reusable Python 3 package to perform **Alanine Scanning Mutagenesis (ASM)** on protein-ligand complexes. It selects binding pocket residues surrounding a ligand based on a distance cutoff, performs virtual alanine mutations, optimizes the mutated structures, and reports the individual binding energy contribution ($\Delta\Delta G$) of each residue along with structural stability changes ($\Delta\text{DOPE}$).

The package leverages **MODELLER** for structure refinement, **PyMOL** for structural manipulations, **OpenBabel** for coordinate conversion, and the **AutoDock Vina** Python API for binding affinity calculations.

---

## Installation

Because the package depends on scientific libraries that require compiled binaries, we recommend using **Conda** to manage dependencies.

### 1. Create a Conda Environment and Install Dependencies

You can create a dedicated environment named `abscan` and install all necessary packages from the `salilab` and `conda-forge` channels.

*Note: MODELLER requires an academic license key. Set the `KEY_MODELLER` environment variable before creating the environment.*

```bash
# Set your MODELLER academic license key (e.g., XXXXXX)
export KEY_MODELLER=XXXXXX

# Create env and install dependencies
conda create -n abscan -c salilab -c conda-forge --yes \
  python=3.12 \
  modeller \
  pymol-open-source \
  vina \
  openbabel \
  matplotlib \
  pandas \
  seaborn \
  biopython
```

### 2. Install the `abscan` Package

Activate the environment and install the package from source:

```bash
# Activate the environment
conda activate abscan

# Install the package in editable mode
pip install -e .
```

---

## Command-Line Usage

Once installed, the pipeline is available via the command-line utility `abscan`.

```bash
abscan -f <protein_ligand.pdb> -n <ligand_residue_number> -d <cutoff_distance> -o <output_directory>
```

### Arguments

| Argument | Long Option | Required | Default | Description |
| :--- | :--- | :---: | :---: | :--- |
| `-f` | `--pdb` | Yes | - | PDB file containing the protein-ligand complex. |
| `-n` | `--resno` | Yes | - | Residue ID (number) of the HETATM ligand in the PDB. |
| `-d` | `--dist` | No | `4.5` | Distance cutoff (Å) from any ligand atom to select pocket residues. |
| `-o` | `--outdir` | Yes | - | Output directory to write all structures and results. |
| `-s` | `--sf` | No | `vina` | Scoring function to use (`vina`, `ad4`, or `vinardo`). |

### Example

Run the scanner on the sample complex `1a4g_A.pdb` (ligand ZMR is residue `466`):

```bash
abscan -f 1a4g_A.pdb -n 466 -d 4.5 -o ./test_run
```

---

## Output Files

The utility will create the output directory and generate the following results:

1. **`abscan_results.csv`**: A CSV file containing:
   - `Residue`: Residue name and original index (e.g. `ARG115`).
   - `DOPE`: Normalized DOPE score of the mutant structure.
   - `dDOPE`: Change in structural stability ($\text{DOPE}_{mutant} - \text{DOPE}_{WT}$). Positive values indicate destabilizing mutations.
   - `Affinity`: Binding affinity calculated by AutoDock Vina (kcal/mol).
   - `ddG`: Change in binding affinity ($\Delta\Delta G_{binding} = \Delta G_{mutant} - \Delta G_{WT}$). Positive values mean the WT residue is stabilizing (contributes positively to binding affinity).
2. **`binding_affinity_changes.png`**: A high-resolution bar plot showing residue contribution to binding affinity ($\Delta\Delta G$).
   - Green bars represent stabilizing residues (affinity decreases when mutated to alanine).
   - Red bars represent destabilizing residues.
3. **`dope_score_changes.png`**: A high-resolution bar plot showing the impact of mutations on protein stability ($\Delta\text{DOPE}$).
   - Red bars indicate destabilizing mutations (increased DOPE).
   - Green bars indicate stabilizing mutations.
4. **Mutant Structures**: All mutated PDB structures (e.g., `ARG_115_ALA.pdb`) and their prepared PDBQT files.

---

## Python API Usage

You can also run the scanner directly in your Python scripts:

```python
from abscan import AlanineScanner, plot_results

# Initialize scanner
scanner = AlanineScanner(
    pdb_file="1a4g_A.pdb",
    resno="466",
    distance=4.5,
    outdir="./test_run",
    scoring_function="vina"
)

# Run the mutagenesis and rescoring workflow
csv_path = scanner.run()

# Generate publication-quality plots
plot_results(csv_path, "./test_run")
```
