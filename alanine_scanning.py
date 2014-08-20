#!/usr/bin/python

import sys
import os
import __main__
__main__.pymol_argv = ['pymol','-rqkc'] # Pymol: quiet and no GUI

#Checking whether the required modules are available in python

#Checking for AutoDockTools installation
try :
	import AutoDockTools
	del AutoDockTools
except ImportError :
	print "AutoDockTools plugin is not available or installed correctly. Please install it from following link:"
	print "http://mgltools.scripps.edu/downloads"
	sys.exit(1)
	pass

#Checking for modeller installtion
try :
	import modeller
	del modeller
except ImportError :
	print "Modeller has not been installed on this system. Please install it from the following link:"
	print "http://salilab.org/modeller/download_installation.html"
	print "If you have already installed modeller, there coule be a problem linking the libraries to python."
	print "This is known to happen on a 64-bit machine. Follow the instructions given on this site:"
	print "http://salilab.org/modeller/9.12/release.html#rpm"
	sys.exit(1)
	pass

#Checking for pymol installation on this system
try :
	import pymol
	del pymol
except ImportError :
	print "Pymol has not been installed on this system. Please install it from the following link:"
	print "http://www.pymol.org/"
	sys.exit(1)
	pass

#Ensuring whether the right number of arguments are being passed with the script
if len(sys.argv) != 5:
	print 'Usage: alanine_scanning.py <protein-ligand complex pdb file> <ligand_resno> <pocket_distance_cutoff> <directory to store results>'
	print 'Ex: alanine_scanning.py 1a4g_A.pdb 466 4.5 test'
	print 'In the above example :'
	print '1a4g_A.pdb is the PDB file containing the protein-ligand complex.'
	print '466 is the residue ID for the ligand - ZMR'
	print '4.5 is the distance cut-off used to select the binding site residues from mentioned ligand atom'
	print 'test is the directory that would be created to store the results.'
	sys.exit(1)

#Ensuring that PDB file exists at the right place 
from os import path, access, R_OK

PDB_FILE = sys.argv[1]
if path.exists(PDB_FILE) and path.isfile(PDB_FILE) and access(PDB_FILE, R_OK):
	pass
else:
	print "Your PDB file does not exist. Please check the path of your PDB file."
	sys.exit(1)

#Checking if the results directory exists
dirname = sys.argv[4]

try:
	os.makedirs(dirname)
except OSError:
	if os.path.exists(dirname):
		print "The directory you specified already exists. Will be using the same directory to write results."
		pass
	else:
		# There was an error on creation, so make sure we know about it
		raise

#Creating a receptor file with only protein atoms with modeller
from os import path, access, R_OK
PDB_FILE = sys.argv[1]

#Checking for the existence of the PDB file and the read permissions
if path.exists(PDB_FILE) and path.isfile(PDB_FILE) and access(PDB_FILE,R_OK):
   pass
else:
   print 'Please check the PDB File entered, either its not accessible or readable'


from modeller import *
from modeller.scripts import complete_pdb
from modeller.optimizers import molecular_dynamics, conjugate_gradients
from modeller.automodel import autosched
log.none()

env = environ()
env.libs.topology.read('${LIB}/top_heav.lib')
env.libs.parameters.read('${LIB}/par.lib')
env.edat.dynamic_sphere=False
env.edat.dynamic_lennard=True
env.edat.contact_shell = 4.0
env.edat.update_dynamic = 0.39

m = complete_pdb(env, sys.argv[1])
m.write(file=sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_receptor.pdb'))

f = open(sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_receptor.pdb'))
lines = f.readlines()
f.close()

f = open(sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_noendreceptor.pdb'),"w")

for line in lines:
  if line!="END"+"\n":
    f.write(line)
f.close()

#Extracting the ligand using pymol from the PDB
import sys, time, os
import pymol

from pymol import cmd

print 'Extracting out the binding site residues with the specified cut-off ='+sys.argv[3]+' Angstroms'
pymol.finish_launching()
pymol.cmd.load(sys.argv[1])
var = ' not'+'(alt \'\'+A)'
pymol.cmd.remove(var)
pymol.cmd.remove('hydrogens')
pymol.cmd.select('ligand', 'resid '+sys.argv[2]+' and hetatm')
pymol.cmd.save(sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_ligand.pdb'), (('ligand')))
pymol.cmd.select('native_pocket', '(byres (resi '+sys.argv[2]+' and hetatm around '+sys.argv[3]+'))')
pymol.cmd.save(sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_native_pocket.pdb'), (('native_pocket')))

native_res_dict = { 'native_res_list' : [] }
pymol.cmd.iterate("native_pocket and name CA and not hetatm","native_res_list.append((resi,resn))",space=native_res_dict)

pymol.cmd.delete('ligand')
pymol.cmd.delete('native_pocket')
pymol.cmd.delete(sys.argv[1].replace('.pdb',''))

#Concatinating the files to get complex.pdb
filenames = [sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_noendreceptor.pdb'),sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_ligand.pdb')]
with open(sys.argv[4]+'/complex.pdb', 'w') as outfile:
    for fname in filenames:
        with open(fname) as infile:
            outfile.write(infile.read())

ligand_dict = { 'lig_list' : [] }

res_dict = { 'res_list' : [] }

#Extracting out binding site residues from pymol
pymol.cmd.load(sys.argv[4]+'/complex.pdb')

pymol.cmd.select('pocket', '(byres (resi '+sys.argv[2]+' and hetatm around '+sys.argv[3]+'))')
pymol.cmd.iterate("pocket and name CA and not hetatm","res_list.append((resi,resn))",space=res_dict)
pymol.cmd.select('ligand', 'resid '+sys.argv[2]+' and hetatm')
pymol.cmd.iterate("ligand","lig_list.append((resn))",space=ligand_dict)
pymol.cmd.save(sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_binding_site.pdb'), (('pocket'))) #Storing the binding site in PDB format
pymol.cmd.quit()

uniq_ligand = []
#print ligand_dict['lig_list']


#Warning if other hetatms are found in the binding site
import re # Importing regular expression
bs_file_open = open(sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_native_pocket.pdb'), "r")
hetatmsite = []

def get_word(text, position):
    words = text.split()
    characters = -1
    for word in words:
        characters += len(word)
        if characters >= position:
            return word

for line in bs_file_open:
	if re.match("(.*)^HETATM(.*)", line):
		hetatm = get_word(line, 12)
		if (hetatm not in hetatmsite):
                   hetatmsite.append(hetatm)

for i in range(0,len(hetatmsite)):
    print "Warning!! HETATM residue - "+ (hetatmsite[i]) + " found at the binding site will not be considered for evaluation of energetics"

for x in ligand_dict['lig_list']:
	if x not in uniq_ligand:
		uniq_ligand.append(x)

#print uniq_ligand[0]

if (len(uniq_ligand) != 1):
	print "Something is wrong with the ligand. Either the resid does not correspond to a ligand or there is more than one ligand"
	sys.exit(1)



binding_site_size = len(res_dict['res_list'])


if (len(res_dict['res_list']) < 1):
		print "No binding site residues were found with this cut-off"
		sys.exit(1)

print 'Running Modeller to create alanine mutants in binding site'

env.io.hetatm = True

def optimize(atmsel, sched):
	#conjugate gradient
	for step in sched:
		step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001) #md 
		refine(atmsel)
		cg = conjugate_gradients()
		cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)

#Performing molecular dynamics to minimize the model
def refine(atmsel):
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = molecular_dynamics(cap_atom_shift=0.39, md_time_step=4.0,
                            md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                (200, 600,
                                 (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                         max_iterations=its, equilibrate=equil)
            init_vel = False

#use homologs and dihedral library for dihedral angle restraints
def make_restraints(mdl1, aln):
   rsr = mdl1.restraints
   rsr.clear()
   s = selection(mdl1)
   for typ in ('stereo', 'phi-psi_binormal'):
       rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
   for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
       rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                spline_dx=0.3, spline_min_points = 5, aln=aln,
                spline_on_site=True)

def alanine_scanning_pdb(pdb_file, residue_id, outfile):
    # Read the original PDB file and copy its sequence to the alignment array:
    mdl1 = model(env, file=pdb_file)
    ali = alignment(env)
    ali.append_model(mdl1, atom_files=pdb_file, align_codes=pdb_file )
   
    #set up the mutate residue selection segment
    s = selection(mdl1.residues[(int(residue_id) - 1)])

    #perform the mutate to alanine operation
    s.mutate(residue_type='ALA')
    
    #get two copies of the sequence.  A modeller trick to get things set up
    ali.append_model(mdl1, align_codes=pdb_file)

    # Generate molecular topology for mutant
    mdl1.clear_topology()
    mdl1.generate_topology(ali[-1])

    # Transfer all the coordinates you can from the template native structure
    # to the mutant (this works even if the order of atoms in the native PDB
    # file is not standard):
    #here we are generating the model by reading the template coordinates
    mdl1.transfer_xyz(ali)

    # Build the remaining unknown coordinates
    mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
    
    #yes model2 is the same file as model1.  It's a modeller trick.
    mdl2 = model(env, file=pdb_file)
   
    #required to do a transfer_res_numb
    #ali.append_model(mdl2, atom_files=modelname, align_codes=modelname)
    #transfers from "model 2" to "model 1"
    mdl1.res_num_from(mdl2,ali)

    #It is usually necessary to write the mutated sequence out and read it in
    #before proceeding, because not all sequence related information about MODEL
    #is changed by this command (e.g., internal coordinates, charges, and atom
    #types and radii are not updated).
    mdl1.write(file=pdb_file+'.tmp')
    mdl1.read(file=pdb_file+'.tmp')

    #set up restraints before computing energy
    #we do this a second time because the model has been written out and read in,
    #clearing the previously set restraints
    make_restraints(mdl1, ali)

    #a non-bonded pair has to have at least as many selected atoms
    mdl1.env.edat.nonbonded_sel_atoms=1

    sched = autosched.loop.make_for_model(mdl1)

    #only optimize the selected residue (in first pass, just atoms in selected
    #residue, in second pass, include nonbonded neighboring atoms)
    #set up the mutate residue selection segment
    s = selection(mdl1.residues[int(residue_id) - 1])

    mdl1.restraints.unpick_all()
    mdl1.restraints.pick(s)

    s.energy()

    s.randomize_xyz(deviation=4.0)

    mdl1.env.edat.nonbonded_sel_atoms=2
    optimize(s, sched)

    #feels environment (energy computed on pairs that have at least one member
    #in the selected)
    mdl1.env.edat.nonbonded_sel_atoms=1
    optimize(s, sched)

    s.energy()

    #give a proper name
    mdl1.write(file=outfile)
    
    #delete the temporary file
    os.remove(pdb_file+'.tmp')

#Now performing alanine scanning for all the residues in the binding site
for i in range(0, binding_site_size):
    alanine_scanning_pdb(sys.argv[4]+'/complex.pdb', int(res_dict['res_list'][i][0]), sys.argv[4]+'/'+str(native_res_dict['native_res_list'][i][1])+'_'+str(native_res_dict['native_res_list'][i][0])+'_ALA_complex.pdb')
    #Creating a receptor file with only protein atoms
    import re # Importing regular expression
    pdb_file_open = open(sys.argv[4]+'/'+str(native_res_dict['native_res_list'][i][1])+'_'+str(native_res_dict['native_res_list'][i][0])+'_ALA_complex.pdb', "r")

    f = open(sys.argv[4]+'/'+str(native_res_dict['native_res_list'][i][1])+'_'+str(native_res_dict['native_res_list'][i][0])+'_ALA.pdb','w')

    for line in pdb_file_open:
         if re.match("(.*)^ATOM(.*)", line):
	    f.write(line)

    f.close()
    pdb_file_open.close()

    print str(native_res_dict['native_res_list'][i][0])+str(native_res_dict['native_res_list'][i][1])+' mutated to ALA.'

print 'Done with generation of mutants.'

#Starting with preparation of the AutodockTools setup
from MolKit import Read
import MolKit.molecule
import MolKit.protein
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation

def receptor_preparation(receptor, output):
	receptor_filename =  receptor
	verbose = None
	repairs = 'checkhydrogens'
	charges_to_add = 'gasteiger'
	preserve_charge_types=None
	cleanup  = "nphs_lps_waters_nonstdres"
	mode = 'automatic'
	delete_single_nonstd_residues = None
	dictionary = None

	mols = Read(receptor)
	mol = mols[0]
	RPO = AD4ReceptorPreparation(mol, mode, repairs, charges_to_add, cleanup, outputfilename=output,
				     delete_single_nonstd_residues=delete_single_nonstd_residues,dict=dictionary)

print "Preparing the Receptor Files ..."
for i in range(0, binding_site_size):
	print 'Preparing the PDBQT file for -'+str(native_res_dict['native_res_list'][i][1])+'_'+str(native_res_dict['native_res_list'][i][0])+'_ALA.pdb'
	receptor_preparation(sys.argv[4]+'/'+str(native_res_dict['native_res_list'][i][1])+'_'+str(native_res_dict['native_res_list'][i][0])+'_ALA.pdb', sys.argv[4]+'/'+str(native_res_dict['native_res_list'][i][1])+'_'+str(native_res_dict['native_res_list'][i][0])+'_ALA.pdbqt')

receptor_preparation(sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_receptor.pdb'),sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_receptor.pdbqt'))

print "Done with preparing all the receptor files."

print "Preparing the ligand file ..."

from AutoDockTools.MoleculePreparation import AD4LigandPreparation
def ligand_preparation(ligand, output):
	ligand_filename =  None
	verbose = None
	add_bonds = False
	repairs = ""
	charges_to_add = 'gasteiger'
	preserve_charge_types=''
	cleanup  = "nphs_lps"
	allowed_bonds = "backbone"
	root = 'auto'
	outputfilename = output
	check_for_fragments = False
	bonds_to_inactivate = ""
	inactivate_all_torsions = False
	attach_nonbonded_fragments = False
	attach_singletons = False
	mode = 'automatic'
	dict = None

	mols = Read(ligand)
	mol = mols[0]
	mol.buildBondsByDistance()
	    
	LPO = AD4LigandPreparation(mol, mode, repairs, charges_to_add, 
			          cleanup, allowed_bonds, root, 
			  	  outputfilename=output,
                                  dict=dict, check_for_fragments=check_for_fragments,
				  bonds_to_inactivate=bonds_to_inactivate, 
				  inactivate_all_torsions=inactivate_all_torsions,
				  attach_nonbonded_fragments=attach_nonbonded_fragments,
				  attach_singletons=attach_singletons)

ligand_preparation(sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_ligand.pdb'), sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_ligand.pdbqt'))

print "Done with preparing the ligand file."

print "Starting the calculations of binding energy"

from PyAutoDock.InternalEnergy import InternalEnergy
from PyAutoDock.MolecularSystem import MolecularSystem
from PyAutoDock.AutoDockScorer import AutoDock41Scorer

def check_types(molecule,std_types):
    d = {}
    for a in molecule.allAtoms:
        d[a.autodock_element] = 1
    mol_types = d.keys()
    non_std = []
    for t in mol_types:
        if t not in std_types:
            non_std.append(t)
    return non_std            

def autodock_scoring(receptor, ligand):
	receptorfilename =  receptor
	ligandfilename =  ligand
	write_file_mode = False
	outputfilename = sys.argv[4]+'/Alanine_Scanning_Binding_Energies_Result.csv'
	parameter_library_filename = None
	exclude_torsFreeEnergy = False
	verbose = None
	ad_scorer = AutoDock41Scorer(exclude_torsFreeEnergy=exclude_torsFreeEnergy)
	supported_types = ad_scorer.supported_types
	receptor = Read(receptorfilename)[0]
	receptor.buildBondsByDistance()
	rec_non_std = ""

	non_std_types = check_types(receptor, supported_types)

	#Checking the format of receptor
	if len(non_std_types):
		rec_non_std = non_std_types[0]
		if len(non_std_types)>1:
			for t in non_std_types[1:]:
				rec_non_std = rec_non_std + '_' + t

	ligand = Read(ligandfilename)[0]
	ligand.buildBondsByDistance()
	lig_non_std = ""
	non_std_types = check_types(ligand, supported_types)

	#Checking the format of ligand
	if len(non_std_types):
		lig_non_std = non_std_types[0]
		if len(non_std_types)>1:
			for t in non_std_types[1:]:
				lig_non_std = lig_non_std + '_' + t

	mode = 'a'

	first = not os.path.exists(outputfilename)
	if write_file_mode:
		mode = 'w'
		first = True

	optr = open(outputfilename, mode)
    	
	if first:
		tstr = "Receptor,Ligand,AutoDock4.1Score,estat,hb,vdw,dsolv,tors\n"
		optr.write(tstr)

	#setup the molecular system
	ostr = ""

	if len(lig_non_std):
		ostr = 'ERROR: unable to score ligand "%s" due to presence of non-standard atom type(s): %s\n' %(ligand.name, lig_non_std)
		optr.write(ostr)

	elif len(rec_non_std):
		ostr = 'ERROR: unable to score receptor "%s" due to non-standard atom type(s): %s\n' %(receptor.name, rec_non_std)
		optr.write(ostr)
	
	else: 
		ms = MolecularSystem()
        	ms.add_entities(receptor.allAtoms)
        	ms.add_entities(ligand.allAtoms)
        	ad_scorer.set_molecular_system(ms)
        	#get the scores
        	#score per term
        	estat, hb,vdw,dsolv = ad_scorer.get_score_per_term()
        	torsEnrg = ligand.TORSDOF * ad_scorer.tors_weight
        	score = estat + hb + vdw + dsolv + torsEnrg
        	ostr = "%19s,%3s,%8.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f\n" %(receptor.name, uniq_ligand[0], score, estat, hb, vdw, dsolv, torsEnrg)
        	optr.write(ostr)

	optr.close()


for i in range(0, binding_site_size):
	autodock_scoring(sys.argv[4]+'/'+str(native_res_dict['native_res_list'][i][1])+'_'+str(native_res_dict['native_res_list'][i][0])+'_ALA.pdbqt', sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_ligand.pdbqt'))

autodock_scoring(sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_receptor.pdbqt'),sys.argv[4]+'/'+sys.argv[1].replace('.pdb','_ligand.pdbqt'))

print 'All the calculations completed successfully. The results are present in '+ sys.argv[4]+'/Alanine_Scanning_Binding_Energies_Result.txt'

correspondences = open(sys.argv[4]+'/correspondences.txt',"w")

for i in range(0,binding_site_size):
    corr = native_res_dict['native_res_list'][i][1]+'_'+native_res_dict['native_res_list'][i][0]+','+res_dict['res_list'][i][1]+'_'+res_dict['res_list'][i][0]+'\n'
    correspondences.write(corr)

correspondences.close()
