# Copyright (c) 2022-2026 Trustees of Boston College
#SPDX-License-Identifier: Apache-2.0
# -*- coding: utf-8 -*-
'''Author@ZeyiHuang
This is the code adapted from Dr.Aaron and Dr. Rachel Kelemen to generate a weblogo showing
the log2 (fold enrichment of point mutation to WT) in a selection.

read_VAMP_para: read the csv parameter file and return WT sequence and highlight annotation information

class LogoFeature: define a class to create feature instances. The instances have 3 attributes: start index, end index and color for highlight

highlight_feature_setting: Use the highlight information from function: read_VAMP_para to generate feature instances. It returns
a list contain all the feature instances.
'''

import pandas as pd
import numpy as np
import logomaker
import matplotlib as plt
import csv

def read_VAMP_para(parameter_file_name):
    p0 = {}

    # Read parameter file (.csv)

    rows = []
    with (open(parameter_file_name, 'r', encoding='utf-8') as parameter_file):
        parameter_file_reader = csv.reader(parameter_file)
        for row in parameter_file_reader:
            rows.append(row)
            # Only read rows with a value in the parameter_name column
            if row[0] != '':
                values = []
                for value in row[1:]:
                    if value != '':
                        values.append(value)
                if len(values) == 0:
                    # Because optional parameters can be left blank
                    p0[row[0]] = None
                else:
                    p0[row[0]] = values

        if p0['WT_sequence'] is None:
            raise ValueError('You did not provide WT sequence.')
        #Make sure that the WT sequence is provided
        else:
            WT_sequence = p0['WT_sequence'][0].upper()
            print('WT sequence: ', WT_sequence)
            p0['WT_sequence'][0] = WT_sequence
        #Provide the WT sequence and make sure that it is in upper case. You can double check the sequence correct or not.
        highlight_cord = {}
        #Create a dict to store the highlight information
        p = {}
        #Create a dict to store adjusted parameters
        for sp in p0.keys():
            if p0[sp]:
                if sp == 'start_index':
                    p[sp] = []
                    for value in p0[sp]:
                        adjusted_value = int(value) -1
                        p[sp].append(adjusted_value)
                elif sp == 'end_index':
                    p[sp] = []
                    for value in p0[sp]:
                        adjusted_value = int(value) -1
                        p[sp].append(adjusted_value)
                elif sp == 'export_file_name':
                    p[sp] = p0['export_file_name'][0]
                else:
                    p[sp] = []
                    for value in p0[sp]:
                        p[sp].append(value)
        if p0['Feature_annotation'][0] == 'TRUE':
            print('Weblogo is featured.')

            for start, end, color in zip (p['start_index'], p['end_index'],p['color']):
                highlight_cord.update({(start, end): color})
            #Extract highlight information for N_tile logo
        else:
            print('Weblogo is not featured.')
        #Create list containing tuples to store pair position information
        natural = []
        artificial = []
        if p['pairs_5prime'] is not None and p['pairs_3prime'] is not None:
            for i in range(len(p['pairs_5prime'])):
                natural.append((int(p['pairs_5prime'][i])-1, int(p['pairs_3prime'][i])-1))

             #Read the pair position if you have pair mutation needed to be processed.

        pairs = natural + artificial
        if pairs != []:
            pairs.sort(key = lambda x: x[0])
            p['pairs'] = pairs
        else:
            p['pairs'] = None

    return p, highlight_cord

class LogoFeature():
#Define a class for generating instances: features. These instance contain the information of highlight region in the weblogo.
    __slots__ = ['begin_index', 'end_index', 'color']

    def __init__(self, begin_index=0, end_index=0, color='b'):
        self.begin_index = int(begin_index)
        self.end_index = int(end_index)
        self.color = color

    def __repr__(self):
        return 'LogoFeture({}, {}, {})'.format(self.begin_index,
                                               self.end_index,
                                               repr(self.color))


def highlight_feature_setting(highlight_cord):
    LogoFeatures = []
    for highlight in highlight_cord:
        start, end = highlight
        LogoFeatures.append(LogoFeature(start, end, highlight_cord[highlight]))
    return LogoFeatures

def make_logo_plot(working_folder_name, fname, base_seq, dat, width=14, height=5, features=None):
    '''
    Make sequence logo plot.

    Parameters
    ----------
    fname: str
        Path to write logo file.
    base_seq: str
        Sequence to show in plot
    dat: DataFrame
        4xN matrix with numbers to show for each base at each position.
    width: int
        Width in inches.
    height: int
        Height in inches.
    features: list
        List of LogoFeatures to add to plot.
    '''
    output_fname_path = working_folder_name + '//' + fname + '.svg'
    x_labels = list()
    base_colors = {'A': 'g', 'C': 'b', 'G': 'darkorange', 'U': 'r'}
    for base, index in zip(base_seq, dat.index.to_list()):
        if (index+1) % 5 == 0:
            x_labels.append('{}\n{}'.format(base, index+1))
        else:
            x_labels.append(base)

    # create logo plot
    logo = logomaker.Logo(dat)
    logo.ax.set_ylabel('$\log_2$ fold enrichment')
    logo.ax.set_xticks(dat.index.to_list())
    logo.ax.set_xticklabels(x_labels)
    logo.style_glyphs_below(flip=False)
    for b, t in zip(base_seq, logo.ax.get_xticklabels()):
        t.set_color(base_colors[b])

    # add feature annotations
    if features is not None:
        for f in features:
            logo.highlight_position_range(pmin=f.begin_index, pmax=f.end_index, color=f.color)

    logo.fig.set_size_inches((width, height))
    logo.fig.savefig(output_fname_path, format='svg')





