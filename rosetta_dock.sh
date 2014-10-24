#!/bin/bash
#This script involves my adaptation of ligand docking as per the rosetta protocol.
#source: https://www.rosettacommons.org/manuals/rosetta3_user_guide/app_ligand_docking.html
#troubleshooting: https://www.rosettacommons.org/content/garbled-ligands-dockingprotocol
#My script is meant for the calculation of ΔΔG during alanine scanning, i.e. the energy difference between a WT protein-ligand complex and an ALA-mutant protein-ligand complex.

#This script does NOT perform alanine scanning or the creation of alanine mutant structures.
#I recommend the use of fixbb.linuxgccrelease for alanine scanning.

#PREREQUISITES:
#rosetta 3.4/3.5.
#Install open babel (on ubuntu/mint distros only):
#$ sudo apt-get install obabel.


#For further queries, please contact
# 1. Deepesh Nagarajan
# Email: 1337deepesh@gmail.com 
# 2. Prof. Nagasuma Chandra
# Email: nchandra@biochem.iisc.ernet.in

#standard usage/error message:
if [[ $# -ne 5 ]]
then
  echo "usage: rosetta_dock.sh <protein.pdb> <ligand.pdb> <preminimization runs> <docking runs> <output_folder>"
  exit
fi

#definitions:
IP_protein=$1
IP_ligand=$2
premin_runs=$3
docking_runs=$4
working_address=$5
anchor=`pwd`

#------------------  END-USER: CHANGE THIS ONLY  ------------------#
  rosetta_directory=~/rosetta3.5
  #installation_directory=~/Desktop/Other_humans/Praveen/rosetta_dock
#------------------------------------------------------------------#

#cleanup
rm -r $working_address
mkdir $working_address

#cleanup input files:
cat $IP_protein|grep "^ATOM  " > $working_address/IP_protein.pdb
cat $IP_ligand|grep "^HETATM" > $working_address/IP_ligand.pdb

cd $working_address

#convert ligand.pdb to ligand.mol:
obabel -ipdb IP_ligand.pdb -omol2 -O IP_ligand.mol2

#convert .mol to .params:
ligand_name=`cat IP_ligand.pdb|head -1|cut -c 18-20`
python $rosetta_directory/rosetta_source/src/python/apps/public/molfile_to_params.py -n $ligand_name -p IP_protein IP_ligand.mol2

if [[ $premin_runs -ne 0 ]]
then
  #repack protein receptor:
  $rosetta_directory/rosetta_source/bin/ligand_rpkmin.linuxgccrelease -database $rosetta_directory/rosetta_database -ex1 -ex2 -ex1aro -extrachi_cutoff 1 -no_optH false -flip_HNQ -docking:ligand:old_estat -docking:ligand:soft_rep -nstruct $premin_runs -s IP_protein.pdb
  for i in {1..10}
  do
    name=`echo $i|awk '{printf"IP_protein_%.4d.pdb\n", $1}'`
    score=`cat $name|grep "^pose "|awk '{print $NF}'`
    echo $name $score
  done|sort -nk2|head -1|awk '{print "cat "$1"|grep \"^ATOM  \"> protein_predock.pdb"}'|bash
  cat IP_protein_0001.pdb >> protein_predock.pdb
else
  #DO NOT repack protein receptor, run as is:
  cat IP_protein.pdb >> protein_predock.pdb
  cat IP_protein_0001.pdb >> protein_predock.pdb
fi

#perform docking:
$rosetta_directory/rosetta_source/bin/ligand_dock.linuxgccrelease -in:file:s protein_predock.pdb -extra_res_fa IP_protein.params -nstruct $docking_runs -database $rosetta_directory/rosetta_database -out:pdb -no_optH false #-packing:use_input_sc

exit
#rosetta.gcc aa 1brs 1 -s 1brs -score -dock -dockFA -dock_score_norepack -ex1 -ex2aro_only -find_disulf -paths ../1brs/paths.txt -nstruct 1 > 1brs_dockFA_score.lo





















