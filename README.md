# cosensus_design
Scripts for consensus protein design using pyRosetta

# Usage

## Design scan will explore all identity thresholds and all consensus thresholds for a given alingment

python3.8 design_scan.py --pdb  --design  --db  --n_designs --csv_out 


--pdb = PDB structure of the protein to be designed

--alignment = Alignment file in clustal format. The first sequence of the alignment should perfetly match the structure

--csv_out = output csv file

--n_designs = number of designs for each combination of identity and consensus thresholds

--design = True or False. If True, residues without consensus will be designed by rosetta design function, if False they will be kept as found in the reference structure

## consensus_module.py contains all the individual functions for bespoke consensus design

get_consensus = generates a consensus sequence given an alignment, the structure and the consensus threshold

consensus_design = Given a structure and a consensus sequence, carries out the consensus design. Can be used with or without the design option. If design = True, all residues without consensus will be designed using Rosetta design


