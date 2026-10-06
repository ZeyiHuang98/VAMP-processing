# Copyright (c) 2022-2026 Trustees of Boston College
#SPDX-License-Identifier: Apache-2.0
# -*- coding: utf-8 -*-
import csv
import datetime

'''
Author: @ZeyiHuang
Finish date: 2022-6-18

Function: Create a single site + single pair randomization library of a sequence. Heading and tailing are constant regions added for PCR amplify. 
Please consider the cloning backbone you are gonna use and design the heading and tailing well.
Please modify the pair_index of your sequence accordingly so that the code can recognize the supposed pair in your sequence. 
'''
#This wild type sequence is MaPyl tRNA wt.
wild_type='GCCCGGATGGTGGAATCGGTAGACACAAGGGATTctaAATCCCTCGGCGTTCGCGCTGTGCGGGTTCAAGTCCCGCTCCGGGTA'


'''This is all pairs, including supposed tertiary pairs in EcTrp tRNA. Keys are the indexes of each pair's 5' base, 
and values are the indexes of each pair's 3' base'''
#Example
EcLeu_pair_index = {0:82,1:81,2:80,3:79,4:78,5:77,6:76,7:13,9:25,10:24,11:23,14:58,17:65,18:66,27:43,28:42,29:41,30:40,31:39,44:57,45:56,46:55,47:54,48:53,
59:75,60:74,61:73,62:72,63:71,64:68}

#address of the file where you want to generate the sequences
output_file = rf'D:\{datetime.date.today()}_EcLeu_VAMP.csv'

#Here is the heading and tailing of EcTrp tRNA.
heading1 = 'CTTGTGGAAAGGACGAAACACC'
tailing1 = 'TTTTTTGCTAGCGGATCGACGAGAgcagagagacctttatttt'

canonical_base = ['A','C','G','T']


'''
Function 1: VAMP_sequence_single_nucleotide_mutation_generation. Create base randomization at every single site.
Function 2: VAMP_sequence_single_pair_mutation_generation. Create base pair randomization at every single base pair you supposed.  
'''

def VAMP_sequence_single_nucleotide_mutation_generation(wild_type, head, tail):
    # Here we create a list, holder_1, in which the element i is the i+1 base of the wild type sequence.
    # We set a collecting list, holder_2, to include all the VAMP sequence. The first sequence would be wild type sequence
    holder_1=[wild_type[e] for e in list(range(len(wild_type)))]
    holder_2=[]
    holder_2.append(head+wild_type+tail)
    # The idea is, for every base in wild type sequence, we convert the base to the other three bases and create a new sequence, which will go into holder_2
    # Here are the dictionaries for converting.
    bases = {'A':0, 'C':1, 'G':2, 'T':3}
    convert={'a': 'A', 'c':'C','g':'G','t':'T' }
    convert_list=['a','c','g','t']
    number=[0,1,2,3]
    for i in list(range(len(wild_type))):
        # For every base of wild type
        for base in number:
            if base != bases[wild_type[i]]:
                #if the base (could be A, C, G or T) is different from the base i+1 in wild type sequence
                holder_1[i]=convert[convert_list[base]]
                # base i+1 is mutated to one of the other three bases. Now the holder_1 becomes a list in which the element i is the i+1 base of the mutated sequence.
                # Finally, the heading and tailing are added
                sequence = head+''.join(holder_1).strip() + tail
                # Create the new mutated sequence and append it into holder_2
                holder_2.append(sequence)
                #reset the holder_1's element back to wild type
                holder_1=[wild_type[e] for e in list(range(len(wild_type)))]

    return holder_2

def VAMP_sequence_single_pair_mutation_generation(wild_type, pair_index, head, tail):
    holder_1 = [wild_type[e] for e in list(range(len(wild_type)))]
    holder_2 = []
    #set up base pairs dictionary. You can add non-canonical pair if you want
    base_pairs = [['A','T'],
                  ['T','A'],
                  ['C','G'],
                  ['G','C'],
                  ['G','T'],
                  ['T','G']]
    for i in list(range(len(wild_type))):
        #check if the base i is the 5' base of a pair
        if i in pair_index:
            for base_pair in base_pairs:
                # skip the wild type pair
                if holder_1[i]==base_pair[0] and holder_1[pair_index[i]]==base_pair[1]:
                    pass
                #all non-wild type pairs are used as mutation
                else:
                    holder_1[i] = base_pair[0]
                    holder_1[pair_index[i]] = base_pair[1]
                    sequence = head + ''.join(holder_1).strip() + tail
                    if sequence not in holder_2:
                        holder_2.append(sequence)
                    holder_1 = [wild_type[e] for e in list(range(len(wild_type)))]
    return holder_2

'''Check and remove the repeat sequence if you run Function 1 and 2 on the same sequences.'''
Checking_singlebase = VAMP_sequence_single_nucleotide_mutation_generation(wild_type, heading1, tailing1)
Checking_singlepair= VAMP_sequence_single_pair_mutation_generation(wild_type, EcLeu_pair_index, heading1,tailing1)
Checking_singlebase_revised = []
print ('Single base lib has ', len(Checking_singlebase), ' members')
print ('Single pair lib has ', len(Checking_singlepair), ' members')
for seq in Checking_singlebase:
    if seq not in Checking_singlepair:
        Checking_singlebase_revised.append(seq)
VAMP = Checking_singlebase_revised+Checking_singlepair
VAMP = sorted(VAMP)
with open(output_file,'w+',encoding='utf-8') as output_file:
    for seq in VAMP:
        output_file.write(seq+'\n')
    output_file.close()