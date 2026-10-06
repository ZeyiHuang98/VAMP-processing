# VAMP-processing
This is a code used to process the high-throughput sequencing data and analysis the mutational profiling result of tRNA. 
Specifically, this project is designed to not only perform statistic analysis of mutational profiling library but also plot the result.
## Authors

- Zeyi Huang

Chatterjee Lab, Boston College

- Aaron J. Maurais

Weerapana Lab, Boston College

## Acknowledgement of Original Work

This code builds on VADER, originally developed by Rachel E. Kelemen
and available at https://github.com/chatterjeelab2022/VADER.

If you use this code, please also cite the original VADER paper:

> Jewel, D., Kelemen, R.E., Huang, R.L. et al. Virus-assisted directed evolution of enhanced suppressor tRNAs in mammalian cells. Nat Methods 20, 95–103 (2023).
> https://doi.org/10.1038/s41592-022-01706-w

## Citation

If you use this code, please cite:
>
## Workflow
#1 Mutational profiling library generation

Open RNA_mutational_profiling_sequence_generate.py. You first need to fill the wild type tRNA sequence at line 14 (upper letter). Second, you need to fill in the pair index of tRNA at line 20. To avoid mistake, it is better to annotate the tRNA's base position. The index of all base pairs will be listed in a 5'base:3'base manner (index = base position -1). You can add tertiary base pair if it is needed. Third, you should replace the header sequence and tailer sequence at line 26,27 that you want to add before and after the tRNA sequence (e.g. for PCR). Then run this .py file and you will receive a .csv file containing all library members, which can be used for ordering oligo pool from Twist Bioscience. For data processing, please generate another .csv file (VAMP-Twist.csv) without head and tail sequence.

#2 Filling data processing parameter file

There are two parameter files which need to be filled. 

First is the Main_VAMP.csv. The filing guidance is basically as same as VADER code reported above. Things to be changed: 1) Row 27: Fill the address of VAMP-Twist.csv. 2) Row 37-38: only regions before and after tRNA are considered as the "constant region". Length can be decided by the users but we recommend at least 10 bases to assure the quality of reads. 3) Row 40-42: all bases of tRNA are considered as randomized bases. the randomized pairs can be filled accordingly. 4) Row 60, 61: the output_format and full_seq_format are replaced by "NNNN...", the numebr of "N" is equal to the tRNA length.

Second is Analysis_VAMP.csv. The wild type sequence will be filled at Row 4: WT_sequence in RNA format. Feature_annotation is default to be true for annotating regions of tRNA with color assigned at Row 6-8. Default annotation is only for stem region. You can change it if you prefer other way to highlight the result in the final weblogo. Row 10 is the name of final .svg Weblogo file. Row 11-12 is the position of all base pairs, listed as 5' and 3' bases positions. 

#3 Run the code

Put all the .fastq files, Main_VAMP.csv, Analysis_VAMP.csv and VAMP-Twist.csv under the working folder (see Row 24 in Main_VAMP.csv). Open Main.py, fill the address of Main_VAMP.csv at line 43 (parameter_file_name), and Analysis_VAMP.csv at Row 87 (VAMP_para). Set Row 82 run_plot_results to be False and Row 86 run_VAMP-analysis to be True. Run Main.py. 

#4 Data plotting

All mutational profiling analysis will be saved under a new "N_tile_logo_data" folder under working folder. For each base or base pair, only mutations at that position (including wild type) will be considered and calculated. All mutations' enrichment will be normalized to wild type and plotted in log2 scale. You will see the Weblogo showing the single base mutation profiling result with color annotated. For base pair result, you can open the base_pair_enrichment_normalized_log2.csv, save it as .xlsx file. The result is listed from top to bottom in a order of base pair position (based on its 5' base) , which can be found in the combined_position_pair_percent_VAMP.csv, Column 1. You can then plot the result using stacking bar graph in excel to generaete the plot. 
