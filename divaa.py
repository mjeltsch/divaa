#!/usr/bin/python3
# -*- coding: UTF-8 -*-
#
# Re-implementation of the DIVAA software as described by
# Rodi, D.J. Mandava, S.; Makowski, L.
# DIVAA: Analysis of Amino Acid Diversity in Multiple Aligned Protein Sequences.
# Bioinformatics 2004, 20, 3481–3489, doi:10.1093/bioinformatics/bth432
# 
# The script takes a single alignment file in fasta format
# and puts out a list of numbers (all between 0 and 1),
# that describe the amino acid diversity at each position
# of the alignment. 1 being a perfectly diverse site (each amino acid has
# the chance of 1/20 to be found at this position) and 0.05 being a
# completely conserved site (chance of 1 to be found at this position).

import sys, os
from Bio import AlignIO
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def run(msa_file):

    alignment = AlignIO.read(open(msa_file), 'fasta')
    sequence_length = alignment.get_alignment_length()

    # Amino acid list (in the order shown below)
    # A C D E F G H I K L M N P Q R S T V W Y - X
    #
    # - stands for a gap in the alignment
    # X stands for anything that was not assigned otherwise
    #
    amino_acid_list = []
    for i in range(sequence_length):
        amino_acid_list.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    
    # Iterate over all columns of the alignment
    for i in range(sequence_length):
        for item in alignment[:, i:i+1]:
            if item.seq == 'A': amino_acid_list[i][0] += 1
            elif item.seq == 'C': amino_acid_list[i][1] += 1
            elif item.seq == 'D': amino_acid_list[i][2] += 1
            elif item.seq == 'E': amino_acid_list[i][3] += 1
            elif item.seq == 'F': amino_acid_list[i][4] += 1
            elif item.seq == 'G': amino_acid_list[i][5] += 1
            elif item.seq == 'H': amino_acid_list[i][6] += 1
            elif item.seq == 'I': amino_acid_list[i][7] += 1
            elif item.seq == 'K': amino_acid_list[i][8] += 1
            elif item.seq == 'L': amino_acid_list[i][9] += 1
            elif item.seq == 'M': amino_acid_list[i][10] += 1
            elif item.seq == 'N': amino_acid_list[i][11] += 1
            elif item.seq == 'P': amino_acid_list[i][12] += 1
            elif item.seq == 'Q': amino_acid_list[i][13] += 1
            elif item.seq == 'R': amino_acid_list[i][14] += 1
            elif item.seq == 'S': amino_acid_list[i][15] += 1
            elif item.seq == 'T': amino_acid_list[i][16] += 1
            elif item.seq == 'V': amino_acid_list[i][17] += 1
            elif item.seq == 'W': amino_acid_list[i][18] += 1
            elif item.seq == 'Y': amino_acid_list[i][19] += 1
            elif item.seq == '-': amino_acid_list[i][20] += 1
            else: amino_acid_list[i][21] += 1

    diversity_list = []
    number_of_sequences = len(alignment)
    for i in range(sequence_length):
        sum_of_probabilities = 0
        for j in range(22):
            sum_of_probabilities += (amino_acid_list[i][j]/number_of_sequences)**2
        diversity_list.append(1/(21*sum_of_probabilities))
    
    # Print the result
    for i, item in enumerate(diversity_list):
        print(f"{i}: {item}")
    
    #
    #plt.plot(diversity_list, color = 'red', linewidth=1, linestyle='dashed')
    #plt.xlabel("amino acid position")
    #plt.ylabel("diversity")
    #plt.savefig(msa_file+".svg")

    df = pd.DataFrame(diversity_list).rolling(20, center = True, min_periods = 1).mean()
    sns_plot = sns.lineplot(data = df, palette = "tab10", linewidth = 1)
    sns_plot.set(xlabel = 'aa position', ylabel = 'amino acid diversity')
    #sns_plot.savefig('figure.png', transparet = True)
    sns_plot.figure.savefig('figure.pdf', dpi = 300)
    #sns_plot.savefig('figure.eps', orientation = 'landscape', dpi = 300)
    sns_plot.figure.savefig('figure.svg')

if __name__ == '__main__':
    msa_file = sys.argv[1]
    if len(sys.argv) == 2 and os.path.isfile(msa_file):
        run(msa_file)
