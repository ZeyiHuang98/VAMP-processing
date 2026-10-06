# Copyright (c) 2022-2026 Trustees of Boston College
#SPDX-License-Identifier: Apache-2.0
# -*- coding: utf-8 -*-

'''@Author: Zeyi Huang
    This is a code to calculate the average enrichment of each mutant at each site/normalized enrichment to WT/log2 normalized enrichment.
    For all the mutants, if the mutants show up in the input but not output, 1 raw count was added to that specific mutant just for later processing
'''

import numpy as np
import math
import os
import csv
import datetime as dt

def log2_fold_enrichment_calculation_N_tile(data_folder_name, Selection, biosample_data_list, WT_sequence):
    '''

    Function to
    Parameters
    ----------
    data_folder_name : the folder path of the final data file stored

    Selection : selection name

    biosample_data_list: a list contain the dicts which store the base fraction information of sequence for each biosample in the selection
                         It contains 3 data set in this order: Input, Output1, Output2
    WT_sequence:  a DNA sequence of wild-type in upper letter. Use it as the reference sequence for normalizing enrichment.

    Returns: return a dict contain log2(normalized fold enrichemnt to WT) for each positon of DNA sequence in the selection
    -------

    '''

    WT_seq =list(WT_sequence)
    bases = ['A','C','G','U']
    position_base_enrichment_ave = {}
    position_base_enrichment_normalized = {}
    position_base_enrichment_normalized_log2  = {}
    missed_base_mutant = []
    for base in bases:
        position_base_enrichment_ave.update({base: [0.0 for _ in range(len(WT_seq))]})
        position_base_enrichment_normalized.update({base: [0.0 for _ in range(len(WT_seq))]})
        position_base_enrichment_normalized_log2.update({base: [0.0 for _ in range(len(WT_seq))]})
    for i in range(len(WT_seq)):
        if biosample_data_list[0][i][WT_seq[i]] == 0:
            raise ValueError('No WT in the pool, no normalization can be done.')
            #if the WT sequence did not show up, then no normalization can be done for VAMP.
        else:
            for base in bases:
                if biosample_data_list[0][i][base] == 0.0:
                    position_base_enrichment_ave[base][i] = 'This mutants does not exist in input.'
                    missed_base_mutant.append(f'{i+1} : {base}')
                else:
                    position_base_enrichment_ave[base][i] = (
                    (biosample_data_list[1][i][base] + biosample_data_list[2][i][base])/(2*biosample_data_list[0][i][base]))

        #Calculating the enrichment of each base at each position
        for base in bases:
            if position_base_enrichment_ave[base][i] == 'This mutants does not exist in input.':
                position_base_enrichment_normalized[base][i] = 1
            else:
                position_base_enrichment_normalized[base][i] = position_base_enrichment_ave[base][i]/position_base_enrichment_ave[WT_seq[i]][i]
        #Normalizing each base fold enrichment to WT bases.
            if position_base_enrichment_normalized[base][i] <= 0:
                position_base_enrichment_normalized_log2[base][i] = 'This number is not logarithmizable'
            else:
                position_base_enrichment_normalized_log2[base][i] = math.log2(position_base_enrichment_normalized[base][i])
        #Logarithmize the normalized fold enrichment

    #Write calculated data into csv for tracing

    ave_data_file = data_folder_name + '//' + f'{dt.date.today()}-{Selection}_position_base_enrichment_ave.csv'
    with open (ave_data_file, 'w', newline = '') as ave:
        keys = list(position_base_enrichment_ave.keys())
        writer = csv.DictWriter(ave, fieldnames = keys)
        writer.writeheader()
        for row in zip(*position_base_enrichment_ave.values()):
            writer.writerow({keys[i]: value for i, value in enumerate(row)})
    normalized_data_file = data_folder_name + '//' + f'{dt.date.today()}-{Selection}_position_base_enrichment_normalized.csv'
    with open (normalized_data_file, 'w', newline = '') as normalized:
        keys = list(position_base_enrichment_normalized.keys())
        writer = csv.DictWriter(normalized, fieldnames = keys)
        writer.writeheader()
        for row in zip(*position_base_enrichment_normalized.values()):
            writer.writerow({keys[i]: value for i, value in enumerate(row)})
    log2_data_file = data_folder_name + '//' + f'{dt.date.today()}-{Selection}_position_base_enrichment_normalized_log2.csv'
    with open (log2_data_file, 'w', newline = '') as normalized_log2:
        keys = list(position_base_enrichment_normalized_log2.keys())
        writer = csv.DictWriter(normalized_log2, fieldnames = keys)
        writer.writeheader()
        for row in zip(*position_base_enrichment_normalized_log2.values()):
            writer.writerow({keys[i]: value for i, value in enumerate(row)})
            

    return position_base_enrichment_normalized_log2, missed_base_mutant

def log2_fold_enrichment_calculation_pair(data_folder_name, Selection, biosample_data_list, WT_sequence, pair_cords):
    """

    Parameters
    ----------
    data_folder_name : the folder path of the final data file stored

    Selection : selection name

    biosample_data_list: a list contain the dicts which store the base fraction information of sequence for each biosample in the selection
                         It contains 3 data set in this order: Input, Output1, Output2
    WT_sequence:  a DNA sequence of wild-type in upper letter. Use it as the reference sequence for normalizing enrichment.

    pair_cords: a list recording the pair position in (5_prime, 3_prime) format

    Returns: return a dict contain log2(normalized fold enrichemnt to WT) for each pair element of DNA sequence (normally RNA element) in the selection
    """
    WT_seq = list(WT_sequence)
    base_pairs = ['AU','UA','GC', 'CG', 'GU', 'UG']
    base_pair_enrichment_ave = {}
    base_pair_enrichment_normalized = {}
    base_pair_enrichment_normalized_log2 = {}
    missed_pair_mutant = []
    for pair in base_pairs:
        base_pair_enrichment_ave.update({pair: [0.0 for _ in range(len(pair_cords))]})
        base_pair_enrichment_normalized.update({pair: [0.0 for _ in range(len(pair_cords))]})
        base_pair_enrichment_normalized_log2.update({pair: [0.0 for _ in range(len(pair_cords))]})
        #create 3 dicts to store the calculated result. Keys are each base pairs, and the values are lists storing
        # the result for each pair element in the sequence
    WT_pairs = []
    for pair_cord in pair_cords:
        p5, p3 = pair_cord[0], pair_cord[1]
        WT_pair = WT_sequence[p5]+WT_sequence[p3]
        WT_pairs.append(WT_pair)
    for i in range (len(pair_cords)):
        if biosample_data_list[0][i][WT_pairs[i]] == 0:
            raise ValueError('No WT in the pool, no normalization can be done.')
            # if the WT sequence did not show up, then no normalization can be done for VAMP.
        else:
            for pair in base_pairs:
                if biosample_data_list[0][i][pair] == 0.0:
                    base_pair_enrichment_ave[pair][i] = 'This mutants does not exist in input.'
                    missed_pair_mutant.append(f'{pair_cords[i][0]+1}-{pair_cords[i][1]+1}: {pair}')
                else:
                    base_pair_enrichment_ave[pair][i] = (
                    (biosample_data_list[1][i][pair] + biosample_data_list[2][i][pair])/(2*biosample_data_list[0][i][pair]))
        # Calculating the enrichment of each pair at each position
        for pair in base_pairs:
            if base_pair_enrichment_ave[pair][i] == 'This mutants does not exist in input.':
                base_pair_enrichment_normalized[pair][i] = 1
            else:
                base_pair_enrichment_normalized[pair][i] = base_pair_enrichment_ave[pair][i]/base_pair_enrichment_ave[WT_pairs[i]][i]
        #Normalizing the enrichment of each pair to the WT pair's enrichment
            if base_pair_enrichment_normalized[pair][i] <= 0:
                base_pair_enrichment_normalized_log2[pair][i] = 'This number is not logarithmizable'
            else:
                base_pair_enrichment_normalized_log2[pair][i] = math.log2(base_pair_enrichment_normalized[pair][i])
        #Logarithmize the normalized fold enrichment
    # Logarithmize the normalized fold enrichment

    # Write calculated data into csv for tracing

    ave_data_file = data_folder_name + '//' + f'{dt.date.today()}-{Selection}_base_pair_enrichment_ave.csv'
    with open(ave_data_file, 'w', newline='') as ave:
        keys = list(base_pair_enrichment_ave.keys())
        writer = csv.DictWriter(ave, fieldnames=keys)
        writer.writeheader()
        for row in zip(*base_pair_enrichment_ave.values()):
            writer.writerow({keys[i]: value for i, value in enumerate(row)})
    normalized_data_file = data_folder_name + '//' + f'{dt.date.today()}-{Selection}_base_pair_enrichment_normalized.csv'
    with open(normalized_data_file, 'w', newline='') as normalized:
        keys = list(base_pair_enrichment_normalized.keys())
        writer = csv.DictWriter(normalized, fieldnames=keys)
        writer.writeheader()
        for row in zip(*base_pair_enrichment_normalized.values()):
            writer.writerow({keys[i]: value for i, value in enumerate(row)})
    log2_data_file = data_folder_name + '//' + f'{dt.date.today()}-{Selection}_base_pair_enrichment_normalized_log2.csv'
    with open(log2_data_file, 'w', newline='') as normalized_log2:
        keys = list(base_pair_enrichment_normalized_log2.keys())
        writer = csv.DictWriter(normalized_log2, fieldnames=keys)
        writer.writeheader()
        for row in zip(*base_pair_enrichment_normalized_log2.values()):
            writer.writerow({keys[i]: value for i, value in enumerate(row)})

    return base_pair_enrichment_normalized_log2, missed_pair_mutant,WT_pairs





