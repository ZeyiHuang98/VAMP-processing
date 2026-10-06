# Copyright (c) 2022-2026 Trustees of Boston College
#SPDX-License-Identifier: Apache-2.0
# -*- coding: utf-8 -*-
'''@Author: Zeyi Huang
    This code calculate a mutagenesis mutation's abundance among all mutants for the corresponding site
    including single base and single pair if needed (like tRNA)
'''
from VAMP_enrichment_calculation import log2_fold_enrichment_calculation_N_tile, log2_fold_enrichment_calculation_pair
import Web_logo
from Frequency3 import base_frequency_cal, base_pair_frequency_cal, pair_graphing_correction, base_graphing_correction
import pandas as pd
import os
import datetime as dt

def VAMP_processing (working_folder_name, Selection, biosamples,start_index,
                     end_index,  date, possible_sequences, VAMP_para_file):
    '''

    Parameters
    ----------
    working_folder_name
    Selection
    biosamples
    start_index
    end_index
    date
    possible_sequences: library sequence list, containing all all possible library sequences
    VAMP_para_file: csv file containing parameter to process mutational profiling data and make plot

    Returns
    -------

    '''
    data_folder_name = working_folder_name + '//N_tile_logo_data'
    p, highlight_cord = Web_logo.read_VAMP_para(VAMP_para_file)
    wild_type = p["WT_sequence"][0]

    pair_base_position = {}
    pairs = {}
    pair_name = []
    if p['pairs'] is not None:
        for pair in p['pairs']:
            x , y = pair
            pair_base_position.update ({str(x):y, str(y):x})
            pairs.update({str(x):y})
            pair_name.append(f'{x+1}-{y+1}')

    if os.path.exists(data_folder_name) == False:
        os.makedirs(data_folder_name)

    combined_base_pair_counts = []

    combined_base_counts = []
    for biosample in biosamples:
        combined_base_counts.append(base_frequency_cal(working_folder_name,biosample, wild_type, start_index, end_index,  date, possible_sequences,pair_base_position))
    combined_base_frequency_cal = base_graphing_correction (wild_type, biosamples, combined_base_counts)
    #Write down the raw data of bases
    with open(rf'{data_folder_name}\{dt.date.today()} combined_position_base_percent_{Selection}.csv', 'w+',
              encoding='utf-8') as output_file:
        output_file.write(f',,{biosamples[0]},{biosamples[1]},{biosamples[2]},'
                          f'{biosamples[0]},{biosamples[1]},{biosamples[2]},'
                          f'{biosamples[0]},{biosamples[1]},{biosamples[2]},'
                          f'{biosamples[0]},{biosamples[1]},{biosamples[2]}\n'
                          f'Base,,A,A,A,C,C,C,G,G,G,U,U,U\n')
        for i in range(len(wild_type)):
            output_file.write(f"{i + 1},{wild_type[i]},{combined_base_frequency_cal[0][i]['A']},"
                              f"{combined_base_frequency_cal[1][i]['A']},"
                              f"{combined_base_frequency_cal[2][i]['A']},"
                              f"{combined_base_frequency_cal[0][i]['C']},"
                              f"{combined_base_frequency_cal[1][i]['C']},"
                              f"{combined_base_frequency_cal[2][i]['C']},"
                              f"{combined_base_frequency_cal[0][i]['G']},"
                              f"{combined_base_frequency_cal[1][i]['G']},"
                              f"{combined_base_frequency_cal[2][i]['G']},"
                              f"{combined_base_frequency_cal[0][i]['U']},"
                              f"{combined_base_frequency_cal[1][i]['U']},"
                              f"{combined_base_frequency_cal[2][i]['U']}\n")  # calculate percentage of each base at each position and output
        output_file.close()
    #generate the data for making N_tile Weblogo
    log2_data, missed_base_mutants = log2_fold_enrichment_calculation_N_tile(data_folder_name, Selection, combined_base_frequency_cal,wild_type)
    df = pd.DataFrame(log2_data)
    LogoFeatures = Web_logo.highlight_feature_setting(highlight_cord)
    Web_logo.make_logo_plot(data_folder_name, p["export_file_name"], p["WT_sequence"][0], df, features=LogoFeatures)

   #If we have pair mutation data need to be processed, make sure that you have fill the pair information in the VAMP_para file
    if p['pairs'] is not None:
        for biosample in biosamples:
            combined_base_pair_counts.append(base_pair_frequency_cal(working_folder_name, biosample,wild_type, start_index, end_index, date,possible_sequences, pairs))
        combined_base_pair_frequency_cal = pair_graphing_correction(wild_type, biosamples,pairs, combined_base_pair_counts)
        log2_pair_data, missed_pair_mutants,WT_pairs = log2_fold_enrichment_calculation_pair(data_folder_name, Selection,combined_base_pair_frequency_cal, wild_type, p['pairs'])

        with open(rf'{data_folder_name}\{dt.date.today()} combined_position_pair_percent_{Selection}_VAMP.csv',
                  'w+', encoding='utf-8') as output_file:
            output_file.write(f',,{biosamples[0]},{biosamples[1]},{biosamples[2]},'
                              f'{biosamples[0]},{biosamples[1]},{biosamples[2]},'
                              f'{biosamples[0]},{biosamples[1]},{biosamples[2]},'
                              f'{biosamples[0]},{biosamples[1]},{biosamples[2]},'
                              f'{biosamples[0]},{biosamples[1]},{biosamples[2]},'
                              f'{biosamples[0]},{biosamples[1]},{biosamples[2]}\n'
                              f'Pair_position,WT_pair,AU,AU,AU,UA,UA,UA,CG,CG,CG,GC,GC,GC,GU,GU,GU,UG,UG,UG\n')
            for i in range(len(pairs)):
                print(pair_name[i])
                output_file.write(f"{pair_name[i]},{WT_pairs[i]},{combined_base_pair_frequency_cal[0][i]['AU']},"
                                  f"{combined_base_pair_frequency_cal[1][i]['AU']},"
                                  f"{combined_base_pair_frequency_cal[2][i]['AU']},"
                                  f"{combined_base_pair_frequency_cal[0][i]['UA']},"
                                  f"{combined_base_pair_frequency_cal[1][i]['UA']},"
                                  f"{combined_base_pair_frequency_cal[2][i]['UA']},"
                                  f"{combined_base_pair_frequency_cal[0][i]['CG']},"
                                  f"{combined_base_pair_frequency_cal[1][i]['CG']},"
                                  f"{combined_base_pair_frequency_cal[2][i]['CG']},"
                                  f"{combined_base_pair_frequency_cal[0][i]['GC']},"
                                  f"{combined_base_pair_frequency_cal[1][i]['GC']},"
                                  f"{combined_base_pair_frequency_cal[2][i]['GC']},"
                                  f"{combined_base_pair_frequency_cal[0][i]['GU']},"
                                  f"{combined_base_pair_frequency_cal[1][i]['GU']},"
                                  f"{combined_base_pair_frequency_cal[2][i]['GU']},"
                                  f"{combined_base_pair_frequency_cal[0][i]['UG']},"
                                  f"{combined_base_pair_frequency_cal[1][i]['UG']},"
                                  f"{combined_base_pair_frequency_cal[2][i]['UG']}\n")  # calculate percentage of each base at each position and output

            output_file.close()
    missed_mutants = missed_pair_mutants + missed_base_mutants
    with open(rf'{data_folder_name}\{dt.date.today()}-{Selection}-missed_VAMP_muatnts.csv','w+', encoding='utf-8') as missed:
        missed.write(f'Missing mutants as shown below:\n')
        for missed_mutant in missed_mutants:
            missed.write(f'{missed_mutant}\n')
        missed.close()
    return 0
