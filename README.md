ABS-Scan
========

Alanine binding site scanning mutagenesis for evaluting the contribution of individual residues at the binding site towards small-molecule ligand recognition.


Web-server
=========
We strongly recommend user-friendly webserver available at <a href="http://proline.biochem.iisc.ernet.in/abscan" target="_blank">http://proline.biochem.iisc.ernet.in/abscan</a>

Graphical output provided on webserver:

<table>
<tr>
<td>
<img src="http://proline.biochem.iisc.ernet.in/abscan/ABSCAN_ddG.png" width="325px" />
</td>
<td>
<img src="http://proline.biochem.iisc.ernet.in/abscan/ABSCAN_residuecontrib.png" width="325px"/>
</td>
</tr>
</table>

Example output can be visualized by clicking here - <a href="http://proline.biochem.iisc.ernet.in/abscan/examples" target="_blank">Examples</a>

Command-line usage
==================
<code>
./alanine_scanning.py -h </br>
usage: alanine_scanning.py [-h] [-f PDBFILE] [-n RESNO] [-d DIST] [-o OUTDIR]</br>

Arguments:</br>

  -h, --help  show this help message and exit
  
  -f PDBFILE  PDB file of protein-ligand complex
  
  -n RESNO    residue number of HETATM in protein-ligand complex
  
  -d DIST     distane-cutoff to use for defining the binding site
  
  -o OUTDIR   output directory to store the results

</code></br>

Ex:<code>./alanine_scanning.py -f 1a4g_A.pdb -n 466 -d 4.5 -o ./test </code>
</br>

In the above example :

1a4g_A.pdb -- It is the PDB file containing the protein-ligand complex.

466 -- is the residue ID for the ligand - ZMR

4.5 -- is the distance cut-off used to select the binding site residues from mentioned ligand atom

test -- is the directory that would be created to store the results.

Dependencies
============
Please ensure following are installed on your system:
<ul>
<li><a href="https://salilab.org/modeller/" target="_blank">Modeller</a></li>
<li><a href="http://mgltools.scripps.edu/downloads" target="_blank">MGL autodockTools</a></li>
<li><a href="www.pymol.org/" target="_blank">Pymol</a></li>
</ul>


