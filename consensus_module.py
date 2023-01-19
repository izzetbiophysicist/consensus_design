from pyrosetta import *
from Bio import AlignIO
from statistics import mode
import argparse
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation

pyrosetta.init()


scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

def pack_relax(pose, scorefxn):

    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.RestrictToRepacking())
    
    # Set up a MoveMapFactory
    mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    mmf.all_bb(setting=True)
    mmf.all_bondangles(setting=True)
    mmf.all_bondlengths(setting=True)
    mmf.all_chi(setting=True)
    mmf.all_jumps(setting=True)
    mmf.set_cartesian(setting=True)

    ## Print informations about structure before apply fast relax
    # display_pose = pyrosetta.rosetta.protocols.fold_from_loops.movers.DisplayPoseLabelsMover()
    # display_pose.tasks(tf)
    # display_pose.movemap_factory(mmf)
    # display_pose.apply(pose)

    fr = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1)
    fr.cartesian(True)
    fr.set_task_factory(tf)
    fr.set_movemap_factory(mmf)
    fr.min_type("lbfgs_armijo_nonmonotone")
    fr.apply(pose)
    return pose


def mutate_repack(pose, posi, amino):
     #Select position to mutate
    mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    mut_posi.set_index(posi)
    
    #Select neighbor positions
    nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_focus_selector(mut_posi)
    nbr_selector.set_include_focus_in_subset(True)
    
    not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)

    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

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
    
    return pose


def design(pose, scorefxn, to_design):
    
    # Initiate selector
    mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    
    for posi in to_design:
        mut_posi.append_index(posi)

    # Select neighbor Positions
    nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_focus_selector(mut_posi)
    nbr_selector.set_include_focus_in_subset(True)
    print(pyrosetta.rosetta.core.select.get_residues_from_subset(nbr_selector.apply(pose)))

    # Select no design area
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
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),not_design))

    # Enable design
    aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    aa_to_design.aas_to_keep("ACDEFHIKLMNQRGPSTVWY")
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))

    # Create Packer
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
    packer.task_factory(tf)

    packer.apply(pose)
    
    return pose


def consensus_design(pose, consensus, scorefxn, design_undef):
    
    
    for position in range(len(consensus)):
        
        if consensus[position] != '-':
            
            print("mutating position"+str(position))
            
            #posi =  pose.pdb_info().number(position+1)
            posi =  position+1
            amino = consensus[position]
            
            pose = mutate_repack(pose, posi, amino)
    
    ### Relax
    
    if design_undef == "yes":
        ################### Design residues with no significant consensus
        
        to_design = [position+1 for position in range(len(consensus)) if consensus[position] == '-']
        
        ### Calculate N gap
        n_gap = 0
        for position in range(len(consensus)):
            if consensus[position] == '-':
                n_gap += 1
                
        if n_gap != 0:
            print('designing')
            pose = design(pose, scorefxn, to_design)
            
            
        ### Second relax
            
    pose = pack_relax(pose, scorefxn)

    
    return pose

def get_consensus(alignment, cv_thresh, pose):
    ### Generate consensus sequence

    al = alignment.get_alignment_length()
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
    
    
    return consensus


    
