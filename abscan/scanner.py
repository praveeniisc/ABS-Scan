import os
import sys
import logging
import pandas as pd
import numpy as np

# Suppress PyMOL startup messages
import __main__
__main__.pymol_argv = ['pymol', '-rqkc']

import pymol
from pymol import cmd

from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import molecular_dynamics, conjugate_gradients
from modeller.automodel import autosched

from .utils import (
    get_ligand_coords,
    get_vina_box,
    prepare_receptor,
    prepare_ligand
)

logger = logging.getLogger("abscan.scanner")

# Optimization helper functions for Modeller
def optimize(atmsel, sched):
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
        refine(atmsel)
        cg = conjugate_gradients()
        cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)

def refine(atmsel):
    md = molecular_dynamics(cap_atom_shift=0.39, md_time_step=4.0, md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                (200, 600, (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                        max_iterations=its, equilibrate=equil)
            init_vel = False

def make_restraints(mdl1, aln):
    rsr = mdl1.restraints
    rsr.clear()
    s = selection(mdl1)
    for typ in ('stereo', 'phi-psi_binormal'):
        rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
    for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
        rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                 spline_dx=0.3, spline_min_points=5, aln=aln, spline_on_site=True)


class AlanineScanner:
    def __init__(self, pdb_file, resno, distance=4.5, outdir="./abscan_results", scoring_function="vina"):
        self.pdb_file = os.path.abspath(pdb_file)
        self.resno = str(resno)
        self.distance = float(distance)
        self.outdir = os.path.abspath(outdir)
        self.scoring_function = scoring_function.lower()
        
        # Validate scoring function
        if self.scoring_function not in ["vina", "ad4", "vinardo"]:
            raise ValueError(f"Unsupported scoring function: {self.scoring_function}. Choose from: vina, ad4, vinardo")
            
        os.makedirs(self.outdir, exist_ok=True)
        
        # Initialize Modeller environment
        log.none()
        self.env = environ()
        self.env.libs.topology.read('${LIB}/top_heav.lib')
        self.env.libs.parameters.read('${LIB}/par.lib')
        self.env.edat.dynamic_sphere = False
        self.env.edat.dynamic_lennard = True
        self.env.edat.contact_shell = 4.0
        self.env.edat.update_dynamic = 0.39
        self.env.io.hetatm = True

    def extract_structure_info(self):
        """Use PyMOL to extract ligand, pocket residues, and complex structures."""
        logger.info("Extracting structures using PyMOL...")
        
        pdb_basename = os.path.basename(self.pdb_file)
        self.ligand_pdb = os.path.join(self.outdir, pdb_basename.replace('.pdb', '_ligand.pdb'))
        self.native_pocket_pdb = os.path.join(self.outdir, pdb_basename.replace('.pdb', '_native_pocket.pdb'))
        self.complex_pdb = os.path.join(self.outdir, "complex.pdb")
        
        # PyMOL operations
        cmd.reinitialize()
        cmd.load(self.pdb_file, "protein_complex")
        
        # Remove alternate conformations (keep only 'A' or empty alt locs)
        cmd.remove("not (alt ''+A)")
        cmd.remove("hydrogens")
        
        # Select ligand
        cmd.select("ligand", f"resid {self.resno} and hetatm")
        ligand_count = cmd.count_atoms("ligand")
        if ligand_count == 0:
            raise ValueError(f"No HETATM residue with number {self.resno} found in PDB file.")
            
        cmd.save(self.ligand_pdb, "ligand")
        
        # Select pocket (residues within distance cutoff, excluding GLY and PRO)
        cmd.select("native_pocket", f"(byres (resi {self.resno} and hetatm around {self.distance} and not resn GLY and not resn PRO))")
        pocket_count = cmd.count_atoms("native_pocket")
        if pocket_count == 0:
            raise ValueError(f"No pocket residues found within {self.distance} Å of ligand (excluding GLY/PRO).")
            
        cmd.save(self.native_pocket_pdb, "native_pocket")
        
        # Get list of native residues (original names/IDs)
        native_res_dict = {'list': []}
        cmd.iterate("native_pocket and name CA and not hetatm", "list.append((resi, resn))", space=native_res_dict)
        self.native_residues = native_res_dict['list']
        
        # Warning if HETATM residues are present in the pocket
        het_dict = {'list': []}
        cmd.iterate("native_pocket and hetatm", "list.append(resn)", space=het_dict)
        uniq_hets = set(het_dict['list'])
        for het in uniq_hets:
            logger.warning(f"HETATM residue '{het}' found in pocket will be excluded from scanning.")
            
        # Clean PyMOL
        cmd.reinitialize()
        
        # Prepare receptor coordinates (Modeller renumbered receptor structure)
        logger.info("Assessing wild-type DOPE score...")
        m = complete_pdb(self.env, self.pdb_file)
        self.wt_dope = m.assess_normalized_dope()
        
        wt_receptor_full_pdb = os.path.join(self.outdir, pdb_basename.replace('.pdb', '_receptor_full.pdb'))
        m.write(file=wt_receptor_full_pdb)
        
        wt_receptor_pdb = os.path.join(self.outdir, pdb_basename.replace('.pdb', '_receptor.pdb'))
        with open(wt_receptor_full_pdb, "r") as infile, open(wt_receptor_pdb, "w") as outfile:
            for line in infile:
                if line.startswith("ATOM  "):
                    outfile.write(line)
        
        # Strip END record from receptor structure for combination
        noend_receptor_pdb = os.path.join(self.outdir, pdb_basename.replace('.pdb', '_noendreceptor.pdb'))
        with open(wt_receptor_pdb, "r") as infile, open(noend_receptor_pdb, "w") as outfile:
            for line in infile:
                if not line.startswith("END"):
                    outfile.write(line)
                    
        # Concatenate protein receptor and ligand to create complex.pdb
        with open(self.complex_pdb, "w") as outfile:
            for fname in [noend_receptor_pdb, self.ligand_pdb]:
                with open(fname, "r") as infile:
                    outfile.write(infile.read())
                    
        # Load complex into PyMOL to map renumbered residues to original ones
        cmd.load(self.complex_pdb, "complex")
        cmd.select("pocket", f"(byres (resi {self.resno} and hetatm around {self.distance} and not resn GLY and not resn PRO))")
        
        renum_res_dict = {'list': []}
        cmd.iterate("pocket and name CA and not hetatm", "list.append((resi, resn))", space=renum_res_dict)
        self.renum_residues = renum_res_dict['list']
        cmd.reinitialize()
        
        if len(self.native_residues) != len(self.renum_residues):
            raise RuntimeError("Mismatch in residue counts between native and renumbered structures.")
            
        logger.info(f"Identified {len(self.native_residues)} residues in the binding pocket for alanine scanning.")

    def mutate_residue(self, renum_id, original_name, original_id):
        """Mutate a single residue to alanine using Modeller and return the DOPE score."""
        logger.info(f"Mutating {original_name} {original_id} (Renumbered ID: {renum_id}) to ALA...")
        
        mutant_complex_pdb = os.path.join(
            self.outdir, f"{original_name}_{original_id}_ALA_complex.pdb"
        )
        
        # Read the complex PDB file
        mdl1 = model(self.env, file=self.complex_pdb)
        ali = alignment(self.env)
        ali.append_model(mdl1, atom_files=self.complex_pdb, align_codes=self.complex_pdb)
        
        # Perform mutation to ALA (residue index is 0-based in Modeller)
        s = selection(mdl1.residues[int(renum_id) - 1])
        s.mutate(residue_type='ALA')
        
        # Alignment tricks to apply mutation
        ali.append_model(mdl1, align_codes=self.complex_pdb)
        mdl1.clear_topology()
        mdl1.generate_topology(ali[-1])
        mdl1.transfer_xyz(ali)
        mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
        
        # Transfer residue numbering from original
        mdl2 = model(self.env, file=self.complex_pdb)
        mdl1.res_num_from(mdl2, ali)
        
        # Write to temporary file and re-read
        tmp_pdb = mutant_complex_pdb + '.tmp'
        mdl1.write(file=tmp_pdb)
        mdl1.read(file=tmp_pdb)
        
        # Restraints and optimization
        make_restraints(mdl1, ali)
        mdl1.env.edat.nonbonded_sel_atoms = 1
        sched = autosched.loop.make_for_model(mdl1)
        
        s = selection(mdl1.residues[int(renum_id) - 1])
        mdl1.restraints.unpick_all()
        mdl1.restraints.pick(s)
        
        s.energy()
        s.randomize_xyz(deviation=4.0)
        
        # Multi-stage optimization
        mdl1.env.edat.nonbonded_sel_atoms = 2
        optimize(s, sched)
        
        mdl1.env.edat.nonbonded_sel_atoms = 1
        optimize(s, sched)
        
        s.energy()
        mdl1.write(file=mutant_complex_pdb)
        
        # Evaluate mutant stability
        dope_score = mdl1.assess_normalized_dope()
        
        # Clean up temp file
        if os.path.exists(tmp_pdb):
            os.remove(tmp_pdb)
            
        # Extract protein atoms from mutated complex to generate receptor structure
        mutant_receptor_pdb = os.path.join(
            self.outdir, f"{original_name}_{original_id}_ALA.pdb"
        )
        with open(mutant_complex_pdb, "r") as infile, open(mutant_receptor_pdb, "w") as outfile:
            for line in infile:
                if line.startswith("ATOM  "):
                    outfile.write(line)
                    
        return dope_score

    def score_vina(self, receptor_pdbqt, ligand_pdbqt, center, box_size):
        """Use Python Vina API to calculate the binding affinity."""
        from vina import Vina
        
        v = Vina(sf_name=self.scoring_function, verbosity=0)
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)
        v.compute_vina_maps(center=center, box_size=box_size)
        
        # Score current pose (without docking)
        scores = v.score()
        return scores[0]  # The first element is the binding affinity (kcal/mol)

    def run(self):
        """Execute the entire alanine scanning mutagenesis workflow."""
        self.extract_structure_info()
        
        # 1. Prepare Ligand (common for all scoring)
        logger.info("Preparing ligand PDBQT...")
        pdb_basename = os.path.basename(self.pdb_file)
        ligand_pdbqt = os.path.join(self.outdir, pdb_basename.replace('.pdb', '_ligand.pdbqt'))
        prepare_ligand(self.ligand_pdb, ligand_pdbqt)
        
        # Calculate box parameters from ligand
        coords = get_ligand_coords(self.ligand_pdb)
        center, box_size = get_vina_box(coords)
        logger.info(f"Ligand center: {center}, Grid Box Size: {box_size}")
        
        results = []
        
        # 2. Score Wild-type (WT) Receptor
        logger.info("Preparing and scoring wild-type receptor...")
        wt_receptor_pdb = os.path.join(self.outdir, pdb_basename.replace('.pdb', '_receptor.pdb'))
        wt_receptor_pdbqt = os.path.join(self.outdir, pdb_basename.replace('.pdb', '_receptor.pdbqt'))
        prepare_receptor(wt_receptor_pdb, wt_receptor_pdbqt)
        
        wt_affinity = self.score_vina(wt_receptor_pdbqt, ligand_pdbqt, center, box_size)
        logger.info(f"Wild-Type binding affinity: {wt_affinity:.4f} kcal/mol")
        
        results.append({
            "Residue": "WT",
            "DOPE": self.wt_dope,
            "dDOPE": 0.0,
            "Affinity": wt_affinity,
            "ddG": 0.0
        })
        
        # 3. Process Mutants
        for i, (orig_id, orig_name) in enumerate(self.native_residues):
            renum_id = self.renum_residues[i][0]
            label = f"{orig_name}_{orig_id}"
            
            try:
                # Modeller mutation
                mutant_dope = self.mutate_residue(renum_id, orig_name, orig_id)
                d_dope = mutant_dope - self.wt_dope
                
                # Prepare mutant receptor PDBQT
                mutant_receptor_pdb = os.path.join(self.outdir, f"{label}_ALA.pdb")
                mutant_receptor_pdbqt = os.path.join(self.outdir, f"{label}_ALA.pdbqt")
                prepare_receptor(mutant_receptor_pdb, mutant_receptor_pdbqt)
                
                # Score mutant complex
                mutant_affinity = self.score_vina(mutant_receptor_pdbqt, ligand_pdbqt, center, box_size)
                
                # ddG = mutant_affinity - wt_affinity
                ddg = mutant_affinity - wt_affinity
                
                logger.info(f"Residue {label} mutated to ALA. Affinity: {mutant_affinity:.4f} kcal/mol, ddG: {ddg:.4f} kcal/mol")
                
                results.append({
                    "Residue": f"{orig_name}{orig_id}",
                    "DOPE": mutant_dope,
                    "dDOPE": d_dope,
                    "Affinity": mutant_affinity,
                    "ddG": ddg
                })
            except Exception as e:
                logger.error(f"Failed to scan residue {label}: {e}")
                
        # 4. Save results to CSV
        results_df = pd.DataFrame(results)
        csv_path = os.path.join(self.outdir, "abscan_results.csv")
        results_df.to_csv(csv_path, index=False)
        logger.info(f"Results successfully saved to {csv_path}")
        
        return csv_path
