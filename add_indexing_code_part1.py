# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 14:54:26 2018

@author: Anton
"""
from Bio import SeqIO 
import os
import glob
from Bio import AlignIO
import numpy as np
import csv
import re
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
#change file names in lines 18,19,21
#This code searches target peptide sequences containing STY phosphosites in each protein against the sequences of those proteins and adds positional info about STY sites in the sequence
file1=open('peptide_atlas_sequences.fasta','r') #file containing sequences of target proteins aka target file
file2=open('Phosphorylation_data_PeptideAtlas_p1_nobrackets.csv','r') #file containing protein accessions and peptide sequences aka template file
file2_reader=csv.reader(file2,delimiter=',')
file3 = open('STY_positions_in_peptide_atlas_proteins_part1.csv', "w",newline='') #file containing STY positions within a certain peptide in a protein
file3.write('Protein Accession'+','+'Peptide sequence with phosphosite'+','+'S positions in protein'+','+'T positions in protein'+','+'Y positions in protein'+'\n')           

accessions=[]
sequences=[]

header_f2=next(file2_reader)

for x in file2_reader: #extract INDIVIDUAL peptide sequences from each protein and the accession of that protein
    accession=x[0]
    sequence=x[2]
    if sequence not in sequences:
        accessions.append(accession)
        sequences.append(sequence)
  
array=np.stack((accessions,sequences),axis=1)

for record in SeqIO.parse(file1, "fasta"): #loop through the FASTA sequences file containing sequences of all peptide Atlas proteins
    serine_positions=[]
    for y in array: #for each of the target PEPTIDE sequences in peptide Atlas file  
        protein=y[0]
        seq=y[1]
        if protein in record.name and seq in record.seq: #if the protein accession matches and if the peptide sequence matches the overall protein sequence
            start_pos_of_peptide_seq=(record.seq).find(seq)  #starting position of the peptide sequence in the protein
            start_pos_of_peptide_seq=start_pos_of_peptide_seq+1 #actual starting position of the peptide sequence in the protein (Python counts from 0)
            S_pos_in_prot=[] #positions of target S residues from target peptide in the protein 
            T_pos_in_prot=[] #positions of target T residues from target peptide in the protein
            Y_pos_in_prot=[] #positions of target Y residues from target peptide in the protein
            S_pos_in_pep=[a.start() for a in re.finditer('S', seq)] #all occurences of S in the peptide sequence and their position in the peptide sequence
            T_pos_in_pep=[a.start() for a in re.finditer('T', seq)] #all occurences of T in the peptide sequence and their position in the peptide sequence
            Y_pos_in_pep=[a.start() for a in re.finditer('Y', seq)] #all occurences of Y in the peptide sequence and their position in the peptide sequence
            if S_pos_in_pep != []: #if there is an S residue in the peptide sequence
                for pos in S_pos_in_pep: #for each position of the S residue in the peptide seqeuence add the starting position of the whole peptide sequence in the protein which gives the position of target S in that protein
                    pos=pos+start_pos_of_peptide_seq
                    S_pos_in_prot.append(pos)
            if T_pos_in_pep != []:
                for pos in T_pos_in_pep:
                    pos=pos+start_pos_of_peptide_seq
                    T_pos_in_prot.append(pos)
            if Y_pos_in_pep != []:
                for pos in Y_pos_in_pep:
                    pos=pos+start_pos_of_peptide_seq
                    Y_pos_in_prot.append(pos)
            #convert all STY position lists to strings so that they can be written into a file. Also get rid of brackets and a comma and apostrophe
            S_pos_in_prot=str(S_pos_in_prot)
            S_pos_in_prot=S_pos_in_prot.replace("'",'')
            S_pos_in_prot=S_pos_in_prot.replace('[','')
            S_pos_in_prot=S_pos_in_prot.replace(']','')
            S_pos_in_prot=S_pos_in_prot.replace(' ','')
            S_pos_in_prot=S_pos_in_prot.replace(',',';')
            T_pos_in_prot=str(T_pos_in_prot)
            T_pos_in_prot=T_pos_in_prot.replace("'",'')
            T_pos_in_prot=T_pos_in_prot.replace('[','')
            T_pos_in_prot=T_pos_in_prot.replace(']','')
            T_pos_in_prot=T_pos_in_prot.replace(' ','')
            T_pos_in_prot=T_pos_in_prot.replace(',',';')
            Y_pos_in_prot=str(Y_pos_in_prot)
            Y_pos_in_prot=Y_pos_in_prot.replace("'",'')
            Y_pos_in_prot=Y_pos_in_prot.replace('[','')
            Y_pos_in_prot=Y_pos_in_prot.replace(']','')
            Y_pos_in_prot=Y_pos_in_prot.replace(' ','')
            Y_pos_in_prot=Y_pos_in_prot.replace(',',';')
            print(protein,seq,S_pos_in_prot,T_pos_in_prot,Y_pos_in_prot)
            file3.write(protein+','+seq+','+S_pos_in_prot+','+T_pos_in_prot+','+Y_pos_in_prot+'\n')
file3.close()
file1.close()
file2.close()
