# cosensus_design
Script using pyRosetta for consensus protein design

# Usage

python3.8 consensus.py --pdb--alignment --out --cons_thresh --design

--pdb = PDB structure of the protein to be designed

--alignment = Alignment file in clustal format. The first sequence of the alignment should perfetly match the structure

--out = output structure file

--cons_thresh = conservation threshold. Residues with conservation above the threshold will be considered as consensus. Residues whose alignment columns have no consensus (i.e. the predominant residue is below the conservation threshold) will be designed using pyRosetta's design routine

--Design = True or False. If True, residues without consensus will be designed by rosetta design function, if False they will be kept as found in the reference structure
