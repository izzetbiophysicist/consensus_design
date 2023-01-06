# cosensus_design
Script using pyRosetta for consensus protein design

# Usage

# Design scan will explore all identity thresholds and all consensus thresholds for a given alingment

python3.8 design_scan.py --pdb  --design  --db  --n_designs --csv_out 


--pdb = PDB structure of the protein to be designed

--alignment = Alignment file in clustal format. The first sequence of the alignment should perfetly match the structure

--csv_out = output csv file

--n_designs = number of designs for each combination of identity and consensus thresholds

--design = True or False. If True, residues without consensus will be designed by rosetta design function, if False they will be kept as found in the reference structure

