#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# This script takes multiple multiple fasta files downloaded from
# NCBI (a search for one specific gene for all othologs)
# and identifies the speies, that is completely represented with
# orthologs for all the genes.
#
# It generates a common file with all fasta files and
# individual multiple fasta files for each input file.
# For the VEGFs and the purpose of this analysis, we
# have downloaded from NCBI https://www.ncbi.nlm.nih.gov/gene/XXXX/ortholog,
# XXXX = 7422, 7423, 7424, 2277, 5228 for VEGF-A, -B, -C, -D and PlGF, respectively.
#
# Notably, if an animal phylum lacks any of these VEGF paralogs,
# NONE of the VEGF sequences from this phylum will be retained!
# For the VEGFs, that means that no birds, crocodiles, amphibia, chrondichthyes, etc.
# are included in the analysis!
#
# Usage: ./1_intersect_multiple_multiple_fastafiles.py uniprot-gene_*.fasta
#

import os, re, sys
from Bio import SeqIO

def extract_organism(seq_record, VERBOSE = True):
    # get organism from square brackets (if there are such)
    re_search_result = re.search('\[(.*)\]', str(seq_record))
    if re_search_result:
        #if VERBOSE: print('Protein sequence detected.')
        organism = re_search_result.group(0)[1:-1]
        # get rid of three-word organisms
        organism = ' '.join(organism.split(' ')[:2])
        if VERBOSE: print('Organism (from protein record {0}): {1}.'.format(seq_record.id, organism))
    # Check whether it is a protein sequence or mRNA sequence
    elif seq_record.id[:3] in ['sp|', 'tr|']:
        #if VERBOSE: print('Protein sequence detected.')
        organism = re.search('OS=(.*) OX=', str(seq_record)).group(1)
        if VERBOSE: print('Organism (from protein record {0}): {1}.'.format(seq_record.id, organism))
    # The rest should be mRNA sequences
    else:
        #if VERBOSE: print('mRNA sequence detected.')
        if 'PREDICTED' in str(seq_record):
            new_str = str(seq_record).split('PREDICTED: ')[1]
        else:
            new_str = seq_record.description.split(' ', 1)[1]
            #print('new_str: {0}'.format(new_str))
        organism = ' '.join(new_str.split(' ', 2)[:2])
        if VERBOSE: print('Organism (from mRNA record {0}): {1}.'.format(seq_record.id, organism))
    return organism

def run(list_of_fasta_files):
    blacklist = ['Ailuropoda melanoleuca']
    # Determine directories of script (in order to load & save the data files)
    APPLICATION_PATH = os.path.abspath(os.path.dirname(__file__))
    FASTA_COMMON_OUTFILE = ''
    # Generate new filename for output
    for file in list_of_fasta_files:
        FASTA_COMMON_OUTFILE += os.path.splitext(file)[0]+'_'
    FASTA_COMMON_OUTFILE = FASTA_COMMON_OUTFILE[:-1]+'_intersection.fasta'
    # Make a list of lists that has as many elements as we have input files
    dictionary_of_lists = {}
    for file in list_of_fasta_files:
        dictionary_of_lists[file] = []
        FASTA_INFILE = file
        records = list(SeqIO.parse(FASTA_INFILE, "fasta"))
        print('{0}: {1} proteins'.format(FASTA_INFILE, len(records)))
        for seq_record in records:
            # The following line contains the actual search for the organism
            # FOR PROTEIN FILES FROM UNIPROT.ORG:
            #organism = re.search('OS=(.*) OX=', str(seq_record)).group(1)
            # FOR mRNA FILES FROM NCBI:
            organism = extract_organism(seq_record, VERBOSE = False)
            if organism not in blacklist:
                if organism not in dictionary_of_lists[file]:
                    dictionary_of_lists[file].append(organism)
                    #print('Organism {0} added...'.format(organism))
        print('Total number of organisms in {0}: {1}'.format(file, len(dictionary_of_lists[file])))
        #print(dictionary_of_lists[file])
    # Convert into sets for intersection operation
    new_set = set(dictionary_of_lists[list_of_fasta_files[0]])
    #print(list_of_fasta_files[0])
    #print(dictionary_of_lists[list_of_fasta_files[0]])
    #print(new_set)
    for i in range(len(list_of_fasta_files)-1):
        new_set = new_set.intersection(set(dictionary_of_lists[list_of_fasta_files[i+1]]))
        #print(new_set)
    print('Common organisms among all files ({0}):'.format(len(new_set)))
    #for item in new_set:
    #    print(item)
    list_of_organisms = list(new_set)
    print('list of organisms:')
    i = 1
    for organism in list_of_organisms:
        print('{0}. {1}'.format(i, organism))
        i += 1
    # Write all entries to a common fasta file
    new_common_records = []
    for file in list_of_fasta_files:
        FASTA_INFILE = file
        FASTA_OUTFILE = os.path.splitext(FASTA_INFILE)[0]+'_intersection.fasta'
        records = list(SeqIO.parse(FASTA_INFILE, "fasta"))
        new_records = []
        # This is more complicated because we want to have all organisms
        # in the same order in both protein and mRNA file
        for organism in list_of_organisms:
            for seq_record in records:
                extracted_organism = extract_organism(seq_record, VERBOSE = False)
                if extracted_organism == organism:
                    new_common_records.append(seq_record)
                    new_records.append(seq_record)
        SeqIO.write(new_records, FASTA_OUTFILE, "fasta")
        print('Output written to {0}'.format(FASTA_OUTFILE))
    # In case we want to have all sequences in one file
    #SeqIO.write(new_common_records, FASTA_COMMON_OUTFILE, "fasta")
    #print('Output written to {0}'.format(FASTA_COMMON_OUTFILE))
    print('Common organisms among all files: {0}.'.format(len(new_set)))

if __name__ == '__main__':
    sys.argv.pop(0)
    file_list_without_intersection = []
    for file in sys.argv:
        if not file.endswith('_intersection.fasta'):
            file_list_without_intersection.append(file)
    run(file_list_without_intersection)
