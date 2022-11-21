from pyrosetta import *
from Bio import AlignIO
from statistics import mode
import argparse

pyrosetta.init()


scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

def consensus_design(pose, consensus, scorefxn):
    
    for position in range(len(consensus)):
        
        if consensus[position] != '-':
            
            print("mutating position"+str(position))
            
            posi = position+1
            amino = consensus[position]
            
            #Select Mutate Position
            mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
            mut_posi.set_index(posi)
            #Select Neighbor Position
            nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
            nbr_selector.set_focus_selector(mut_posi)
            nbr_selector.set_include_focus_in_subset(True)
            # Select No Design Area
            not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)
            # The task factory accepts all the task operations
            tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
            # These are pretty standard
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
        
            # Disable Packing
            prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
            prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
            tf.push_back(prevent_subset_repacking)
        
            # Disable design
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),not_design))
        
            # Enable design
            aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
            aa_to_design.aas_to_keep(amino)
            tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))
        
            # Create Packer
            packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
            packer.task_factory(tf) 
            packer.apply(pose)
    
    
        #FastRelax Protocol
    mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    mmf.all_bb(setting=True)
    mmf.all_bondangles(setting=True)
    mmf.all_bondlengths(setting=True)
    mmf.all_chi(setting=True)
    mmf.all_jumps(setting=True)
    mmf.set_cartesian(setting=True)
    
    fr = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1)
    fr.cartesian(True)
    fr.set_task_factory(tf)
    fr.set_movemap_factory(mmf)
    fr.min_type("lbfgs_armijo_nonmonotone")
      
    #### First relax
    fr.apply(pose)
    
    ################### Design residues with no significant consensus
    
    mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    for position in range(len(consensus)):
        if consensus[position] == '-':
            posi = position+1
            mut_posi.append_index(posi)
            #### Print selected residues for mut posi
            #print(pyrosetta.rosetta.core.select.get_residues_from_subset(mut_posi.apply(pose)))

    # Select Neighbor Position
    nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_focus_selector(mut_posi)
    nbr_selector.set_include_focus_in_subset(True)
    #print(pyrosetta.rosetta.core.select.get_residues_from_subset(nbr_selector.apply(pose)))

    # Select No Design Area
    not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)
    #### Print residues to NOT design
    #print(pyrosetta.rosetta.core.select.get_residues_from_subset(not_design.apply(pose)))

    # The task factory accepts all the task operations
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

    # These are pretty standard
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    # Disable Packing
    prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
    prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
    tf.push_back(prevent_subset_repacking)

    # Disable design
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),not_design))

    # Enable design
    aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    aa_to_design.aas_to_keep("ACDEFGHIKLMNPQRSTVWY")
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))

    # Create Packer
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
    packer.task_factory(tf)
    
    packer.apply(pose)
    
    ### Second relax
    
    fr.apply(pose)
        
    return pose

if __name__ == "__main__":
   
    parser = argparse.ArgumentParser(description='')
    
    parser.add_argument('--pdb', type=str, required=True)
    parser.add_argument('--alignment', type=str, required=True)
    parser.add_argument('--out', type=str, required=True)
    parser.add_argument('--cons_thresh', type=str, required=True)

    args = parser.parse_args()

    alignment = AlignIO.read(open(args.alignment), "clustal")

    al = alignment.get_alignment_length()

    #### Conservation threshold
    cv_thresh = float(args.cons_thresh)

    ### Generate consensus sequence
    consensus = []
    for pos in range(al):
        column = []
        if alignment[0][pos] != '-':
            for align in alignment:
                if align[pos] != '-':
                    column.append(align[pos])
            ### verify percentage and thresholds
            mode_count = 0
            for i in range(len(column)):
                if column[i] == mode(column) :
                    mode_count = mode_count + 1
                conservation = mode_count/len(column)
            
            if conservation  > cv_thresh:
                consensus.append(mode(column))
            else:
                consensus.append('-')
       





    pose = pose_from_pdb(args.pdb)
    pose = consensus_design(pose, consensus, scorefxn)    


    pose.dump_pdb(args.out)


    



    


