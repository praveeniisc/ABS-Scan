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
<img src="http://proline.biochem.iisc.ernet.in/abscan/ABSCAN_ddG.png" width="350px" />
</td>
<td>
<img src="http://proline.biochem.iisc.ernet.in/abscan/ABSCAN_residuecontrib.png" width="350px"/>
</td>
</tr>
</table>

Command-line usage
==================
<code>
Usage: alanine_scanning.py &lt;protein-ligand complex pdb file&gt; &lt;ligand_resno&gt; &lt;pocket_distance_cutoff&gt; &lt;directory to store results&gt;
</code></br>

Ex:<code> ./alanine_scanning.py 1a4g_A.pdb 466 4.5 test</code>
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


