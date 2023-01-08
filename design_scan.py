#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 18:30:49 2022

@author: lucas
"""


from consensus_module import *
import pandas as pd
import os
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from pyrosetta import *
from Bio import AlignIO
from Bio import SeqIO



pyrosetta.init()

os.system("mkdir ./structures")

## Params
def consensus_steps(id_thresh, cv_thresh, design, pose):
            
    ### read blast results
    blast = pd.read_table("blastp.txt", header=None)
    overcutoff = [x for x in range(len(blast[2])) if blast[2][x] >= id_thresh if blast[2][x] != 100]
    
    
    ## add filtered sequences to input file
    scount=0
    with open("family.fasta", 'r') as infile, open("overcutoff.fasta", 'w') as outfile:
      for line in infile:
        if line.startswith(">"):
            scount += 1
            if scount in overcutoff:
                outfile.write(line)
        else:
            if scount in overcutoff:
                outfile.write(line)
        
    os.system("cat seedprot.fasta overcutoff.fasta > input.fasta")
    os.system("clustalo --in=input.fasta --out=family.aln --force --outfmt=clustal --wrap=80")
    
    
    ##### Call consensus
    
    
    alignment = AlignIO.read('family.aln', "clustal")
    
    
    al = alignment.get_alignment_length()
    
        
        
    ### Generate consensus sequence
    consensus = get_consensus(alignment, cv_thresh, pose)
        
    pose_consensus = consensus_design(pose, consensus, scorefxn, design)    
    
    
    
    #### create pandas dataframe with the whole scan
    
    design_result = pd.DataFrame({'consensus sequence':''.join(consensus),'final sequence':pose_consensus.sequence(),'consensus threshold':cv_thresh, 'identity threshold':id_thresh, 'number of sequences':len(alignment), 'rosetta score':scorefxn(pose_consensus), 'design':str(design), 'design number':'NA'}, index=[0])
        
    return design_result, pose_consensus

def run_blast(querry, db):
    os.system("makeblastdb -in "+db+" -dbtype prot")
    os.system("blastp -query "+querry+" -db family.fasta -evalue 1e-6 -num_threads 4 -out blastp.txt -outfmt 6")

def write_sequence_from_pose(pose, out):
    seq = '>seedprot\n'+pose.sequence()+'\n'
    f = open(out, 'wb')
    f.write(seq.encode('utf-8'))
    f.close()


if __name__ == "__main__":
   
    parser = argparse.ArgumentParser(description='')
    
    parser.add_argument('--pdb', type=str, required=True)
    parser.add_argument('--db', type=str, required=True)  ### multifasta for building blast database
    parser.add_argument('--design', type=str, required=True) ## design or not
    parser.add_argument('--n_designs', type=int, required=True) ## number of designs to generate for each combinations of parameters
    parser.add_argument('--csv_out', type=str, required=True) ## number of designs to generate for each combinations of parameters


    args = parser.parse_args()

    ## take arguments    
    db = args.db
    design = args.design
    n_designs = args.n_designs
    csv_file = args.csv_out
    
    ### Load pose
    pose = pose_from_pdb(args.pdb)

    
    ### run blast
    querry='seedprot.fasta'

    write_sequence_from_pose(pose, querry)
    run_blast(querry, db)
    
    ## Create dataframe with results

    ### run scan
    cons_vec = [round(x*0.1,1) for x in range(11)]
    id_vec = [x*10 for x in range(11)]                
    final_data = pd.DataFrame({'consensus sequence':'original structure','final sequence':pose.sequence(),'consensus threshold':'original structure', 'identity threshold':'relaxed original structure', 'number of sequences':'relaxed original structure', 'rosetta score':scorefxn(pose), 'design':'original structure', 'design number':'original structure'}, index=[0])

    
    for design_number in range(n_designs):
        
        #### reloading pose each loop because in previous tests the pose object is being modified - find leakage
        
        pose = pose_from_pdb(args.pdb)
        
        pose_relax = pack_relax(pose, scorefxn)
        pose_relax.dump_pdb('./structures/relax_'+str(design_number)+'.pdb')
        #pose_relax = pose_from_pdb('./structures/relax_'+str(n_designs)+'.pdb')

        
        relax_data = pd.DataFrame({'consensus sequence':'relaxed original structure','final sequence':pose_relax.sequence(),'consensus threshold':'relaxed original structure', 'identity threshold':'relaxed original structure', 'number of sequences':'relaxed original structure', 'rosetta score':scorefxn(pose_relax), 'design':'relaxed original structure', 'design number':'relaxed original structure '+ str(design_number)}, index=[0])
        final_data = final_data.append(relax_data)

        for cv_thresh in cons_vec:
            for id_thresh in id_vec:
                    
                ## relax and dump
                
                consensus_result = consensus_steps(id_thresh, cv_thresh, design, pose_relax)
                pose_consensus = consensus_result[1]
                
                pose_consensus.dump_pdb('./structures/design_cons_'+str(cv_thresh)+'_id_'+str(id_thresh)+"_design_"+str(design)+"_experiment_"+str(design_number)+".pdb") ##### NAME POSE
    
                new_entry = consensus_result[0]
                
                new_entry['design number'] = design_number
                
                final_data = final_data.append(new_entry)
                final_data.to_csv('tmp_dataframe.csv')
        
    final_data.to_csv(csv_file)
