# Copyright (c) 2022-2026 Trustees of Boston College
#SPDX-License-Identifier: Apache-2.0
# -*- coding: utf-8 -*-
"""
@author: ZeyiHuang
"""
import numpy as np
import csv


def base_pair_frequency_cal(working_folder_name, biosample, wild_type,start_index, end_index,  date, possible_sequences,pairs):
    filter_file_name = rf'{working_folder_name}\{date} {biosample} mismatch filtered.txt'
    base_pairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']
    tRNA_numbered = list(np.arange(start_index, end_index+1))
    with open(filter_file_name, 'r', encoding='utf-8') as filter_file:
        holder = filter_file.readlines()
        variant_mark = [False for _ in range(len(holder))]
        sequence_count = 0
        Pool_sequence_count = 0
        position_pair_counts = [{'AU': 0, 'UA': 0, 'CG': 0, 'GC': 0, 'GU': 0, 'UG': 0, } for _ in range(len(tRNA_numbered))]
        for sequence in holder:
            rna_read = sequence.upper().replace('T', 'U')
            tRNA_sequence = [rna_read[e] for e in tRNA_numbered]  # extract tRNA from sequence result
            for i in range(len(tRNA_sequence)):
                if str(i) in pairs:
                    if tRNA_sequence[i] + tRNA_sequence[pairs[str(i)]] in base_pairs:
                        if tRNA_sequence[i] + tRNA_sequence[pairs[str(i)]] != wild_type[i] + wild_type[pairs[str(i)]]:
                            if ''.join(tRNA_sequence).strip() in possible_sequences:
                                position_pair_counts[i][tRNA_sequence[i] + tRNA_sequence[pairs[str(i)]]] += 1  # make sure that here using the library only VAMP
                                # pair mutation check. If it is a pair mutation then add the count.
                                Pool_sequence_count +=1
            if ''.join(tRNA_sequence).strip() != wild_type:
                variant_mark[sequence_count] = True
            sequence_count += 1
        wild_type_count = variant_mark.count(False)  # wild type sequence counts
        Pool_sequence_count += wild_type_count
        print('wild_type_count:', wild_type_count)
        print('In base pair library pool, wild_type occupies: ', wild_type_count/Pool_sequence_count)
        for i in range(len(tRNA_numbered)):
            if str(i) in pairs:
                if wild_type[i] + wild_type[pairs[str(i)]] in base_pairs:
                    position_pair_counts[i][wild_type[i] + wild_type[pairs[str(i)]]] += wild_type_count
        filter_file.close()

    return position_pair_counts

def base_frequency_cal(working_folder_name,biosample, wild_type, start_index, end_index,  date, possible_sequences,pair_base_position):
    filter_file_name = rf'{working_folder_name}\{date} {biosample} mismatch filtered.txt'
    tRNA_numbered = list(np.arange(start_index, end_index+1))
    with open(filter_file_name, 'r', encoding='utf-8') as filter_file:
        holder = filter_file.readlines()
        variant_mark = [False for _ in range(len(holder))]
        sequence_count = 0
        Pool_sequence_count = 0
        position_base_counts = [{'A': 0, 'C': 0, 'G': 0, 'U': 0} for _ in range(
            len(tRNA_numbered))]  # set a list containing dicts, each dicts encode every single position base mutation counts
        for sequence in holder:
            rna_read = sequence.upper().replace('T', 'U')
            tRNA_sequence = [rna_read[e] for e in tRNA_numbered]  # extract tRNA from sequence result
            for i in range(len(tRNA_sequence)):
                if str(i) in pair_base_position:
                    if tRNA_sequence[i] != wild_type[i] and tRNA_sequence[pair_base_position[str(i)]] == wild_type[
                        pair_base_position[str(i)]]:
                        if ''.join(tRNA_sequence).strip() in possible_sequences:
                            position_base_counts[i][tRNA_sequence[i]] += 1# distinguish the mutation is on the pair region or in a non-pair region.
                                                                          #If it is a pair, make sure that it is not a double-side point mutation
                                                                          # e.g. could be AU to GU, but not AU to GC
                            Pool_sequence_count += 1
                elif tRNA_sequence[i] != wild_type[i]:  # if not wild type base, then + 1 to this position's this base's counts
                    if ''.join(tRNA_sequence).strip() in possible_sequences:
                        position_base_counts[i][tRNA_sequence[i]] += 1 # Here is the single base region mutation
                        Pool_sequence_count += 1
            if ''.join(tRNA_sequence).strip() != wild_type:
                variant_mark[sequence_count] = True
            sequence_count += 1
        wild_type_count = variant_mark.count(False)  # wild type sequence counts
        Pool_sequence_count += wild_type_count
        print('wild_type_count:', wild_type_count)
        print('In base library pool, wild_type occupies: ', wild_type_count/Pool_sequence_count)
        for i in range(len(tRNA_numbered)):
            position_base_counts[i][wild_type[i]] += wild_type_count
        filter_file.close()



    return position_base_counts

def pair_graphing_correction(wild_type, biosamples,pairs,combined_base_pair_counts):

    base_pairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']

    combined_base_pair_frequencies = []
    for counts in range(len(combined_base_pair_counts)):
        position_pair_percent = [{'AU': 0.0, 'UA': 0.0, 'CG': 0.0, 'GC': 0.0, 'GU': 0.0, 'UG': 0.0} for _ in
                                 range(len(pairs))]
        p = 0
        for i in range(len(wild_type)):
            if str(i) in pairs:
                position_pair_total = (combined_base_pair_counts[counts][i]['AU'] + combined_base_pair_counts[counts][i]['UA'] + combined_base_pair_counts[counts][i]['CG']
                                       + combined_base_pair_counts[counts][i]['GC'] + combined_base_pair_counts[counts][i]['GU'] + combined_base_pair_counts[counts][i]['UG'])
                if position_pair_total == 0:
                    raise ValueError(f'No expected sequence found for pair {i + 1}-{pairs[str(i)] + 1}')
                if biosamples[counts] != biosamples[0]:
                    for base_pair in base_pairs:
                        if combined_base_pair_counts[0][i][base_pair]!= 0 and combined_base_pair_counts[counts][i][base_pair] == 0:
                            combined_base_pair_counts[counts][i][base_pair] = 1
                            position_pair_total += 1
                for base_pair in base_pairs:
                    position_pair_percent[p][base_pair] = combined_base_pair_counts[counts][i][base_pair] / position_pair_total
                p += 1
        combined_base_pair_frequencies.append(position_pair_percent)
    return combined_base_pair_frequencies

def base_graphing_correction(wild_type, biosamples, combined_base_counts):
    bases = ['A', 'C', 'G', 'U']
    combined_base_frequencies = []
    for counts in range(len(combined_base_counts)):

        position_base_percent = [{'A': 0.0, 'C': 0.0, 'G': 0.0, 'U': 0.0, } for _ in range(len(wild_type))]
        for i in range(len(wild_type)):
            position_base_total = combined_base_counts[counts][i]['A'] + combined_base_counts[counts][i]['C'] + combined_base_counts[counts][i]['G'] + combined_base_counts[counts][i]['U']
            if position_base_total == 0:
                raise ValueError(f'No expected sequence found for position {i + 1}')
            if biosamples[counts] != biosamples[0]:
                for base in ['A', 'C', 'G', 'U']:
                    if combined_base_counts[0][i][base] != 0 and combined_base_counts[counts][i][base] == 0 :
                        combined_base_counts[counts][i][base] = 1
                        position_base_total += 1
            for base in bases:
                position_base_percent[i][base] = combined_base_counts[counts][i][base] / position_base_total
        combined_base_frequencies.append(position_base_percent)
    return combined_base_frequencies
