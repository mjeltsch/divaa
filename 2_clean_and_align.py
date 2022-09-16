#!/usr/bin/python3
# -*- coding: UTF-8 -*-
#

import re, os
from datetime import datetime
from shutil import copyfile, rmtree, which
from time import sleep
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqFeature
from phylolib import execute_subprocess
from VEGFA_isoform_list import load_isoform_dict, get_manual_overrides

# This script takes a fasta dump from uniprot.org for one specific gene
# (in this implementation, VEGFC). The uniprot advanced search is used
# for Gene name "VEGFC" (which yielded on May 3rd 413 entries). This list
# is downloaded (format: uncompressed fasta including isoforms), and further
# cleaned by deleting partial sequences and splice isoforms using the cutoffs
# defined for each individual gene (e.g. for VEGFC, sequences are retained
# if they are longer than 375 amino acid residues and shorter than 425 amino
# acid residues). When multiple sequences are available for the same species,
# the shorter sequences are removed as they represent mostly splice isoforms
# or N-terminally incomplete entries (manual examination of all these cases).
# 
# Sequences can also be excluded manually. E.g. one sequence was manually
# excluded for VEGFC due to a likely faulty gene prediction (Mandrillus
# leucophaeus).
#
# From the 65 species in the intersection (which have sequences for all 5 VEGFs)
# a few need to be removed since some sequences are clearly bogus,notably:
#
# VEGF-D Ailuropoda melanoleuca (Giant Panda)
#

def return_organism(rec):
    return rec.organism

def return_length(rec):
    return -len(rec.seq)

def return_id(rec):
    return rec.id

def print_and_log_message(message):
    global logfile
    with open(logfile, 'a+') as appendfile:
        appendfile.write(message+'\n')
    print(message)

def clean_sequences(infile, gene, data):
    # data[1] = minimal length to retain a sequence
    # data[2] = maximal length to retain a sequence
    # data[3] = blacklist

    # Clean the output directory
    if os.path.isdir(output_directory):
        rmtree(output_directory)
    os.mkdir(output_directory)

    records = list(SeqIO.parse(infile, 'fasta'))
    print_and_log_message('Reading in {0} {1} sequences.'.format(len(records), gene))
    tmp_records = []
    list_of_mutiple_sequences = []

    # Remove wrong isoforms from VEGFA
    if gene[0:5] == 'VEGFA':
        isoform_dict = load_isoform_dict()
        manual_isoform_overrides = get_manual_overrides()
        for sequence_id, isoform_data in manual_isoform_overrides.items():
            isoform_dict[sequence_id] = isoform_data
        isoform = gene[5:8]
        for seq_record in records:
            if isoform_dict[seq_record.id][0] == isoform:
                tmp_records.append(seq_record)
            # Take also the hypothetical N-terminally extended L isoforms
            elif isoform_dict[seq_record.id][0] == isoform+'L':
                tmp_records.append(seq_record)
        records = tmp_records
        print_and_log_message('Remaining {0} sequences after removing all isoforms other than {1}:'.format(len(records), isoform))
        print_and_log_message('{0}:'.format(records))

    # Remove too long (only a few) and too short (truncated) sequences
    tmp_records = []
    for seq_record in records:
        # Let them be as long as possible....
        #if len(seq_record.seq) > data[1] and len(seq_record.seq) < data[2]:
        if len(seq_record.seq) > data[1]:
            #print(seq_record.id)
            if seq_record.id not in data[3]:
                tmp_records.append(seq_record)
            else:
                print_and_log_message('{0} is blacklisted and will be removed.'.format(seq_record.id))
    records = tmp_records
    print_and_log_message('Remaining {0} sequences after removing too short (<{1}) and too long (>{2}) and blacklisted sequences: {3}:'.format(gene, data[1], data[2], len(records)))
    print_and_log_message('{0}:'.format(records))

    for seq_record in records:
        result = re.search(' OS=(.+?) OX=', seq_record.description)
        if result:
            seq_record.organism = result.group(1)
        else:
            seq_record.organism = 'Unknown organism'

    # This removes identical sequences (i.e. aa sequence is identical and organism is identical)
    # List of sequence/organism
    record_seqs = list()
    # New list to write to file
    new_records = list()
    for record in records:
        if record.seq not in record_seqs:
            record_seqs.append([record.seq, record.organism])
            new_records.append(record)
    SeqIO.write(new_records, final_fasta_file_without_identicals, 'fasta')
    print_and_log_message('Remaining {0} sequences after removing identical sequences (if aa sequence and organism are the same): {1}:'.format(gene, len(new_records)))
    print_and_log_message('{0}:'.format(new_records))

    # If a species has more than one sequence, this function gathers all of them into
    # one file and aligns them for manual comprison. These files are prefixed with "1_"
    #
    for seq_record in new_records:
        # Extract the scientific ("latin") organism name
        result = re.search(' OS=(.+?) OX=', seq_record.description)
        if result:
            seq_record.organism = result.group(1)
        else:
            seq_record.organism = 'Unknown organism'
        organism_file = output_directory+'/'+seq_record.organism.replace(' ', '_') + '.fasta'
        # Check whether there is already a sequence from this species
        # If yes: - rename the file (add .tmp ending)
        #         - add file path to the list that keeps track of species with multiple sequences
        if os.path.isfile(organism_file):
            copyfile(organism_file, organism_file+'.tmp')
            if organism_file not in list_of_mutiple_sequences:
                list_of_mutiple_sequences.append(organism_file)
        # Write the new sequence to a fasta file with canonical file name
        SeqIO.write(seq_record, organism_file, 'fasta')
        # Append the contents of the renamed (.tmp) file to the end of the new fasta file
        # Somehow the sleep is necessary if file IO is slow....
        sleep(0.5)
        if os.path.isfile(organism_file+'.tmp'):
            with open(organism_file, 'a+') as appendfile:
                with open(organism_file+'.tmp', 'r') as infile:
                    appendfile.write(infile.read())
        # Remove the .tmp file
            os.remove(organism_file+'.tmp')
    # Alignment with muscle to select the longest isoform if multiple isoforms are present in one species
    for organism_file in list_of_mutiple_sequences:
        organism_file.split('/')[1][:-6]
        multiple_fasta_file = '{0}/1_{1}_multiple.fasta'.format(gene, organism_file.split('/')[1][:-6])
        aligned_fasta_file = '{0}/1_{1}_aligned.fasta'.format(gene, organism_file.split('/')[1][:-6])
        os.rename(organism_file, multiple_fasta_file)
        # Align multiple fasta file (silently)
        output, error = execute_subprocess('Running muscle...', 'muscle -in {0} -out {1}'.format(multiple_fasta_file, aligned_fasta_file), verbose=False)
        #print(output, error)

        # Try to identify the best isoform for the alignment
        # Only VEGFA has isoform data at the moment!
        # 
        single_species_records = list(SeqIO.parse(multiple_fasta_file, 'fasta'))
        for record in single_species_records:
            # for VEGFA only
            if gene[0:5] == 'VEGFA':
                isoform = isoform_dict[record.id][0]
                if 'L' in isoform: i = 0
                elif 'B' in isoform: i = 1
                else: i = 2
                # Insert gaps (%)
                record.features.append(isoform_dict[record.id][1]) # numeric
            # all other VEGFs, there is not isoform dictionary
            else:
                # Insert "dummy" data for all VEGFs except VEGFA
                # Isoform = 0
                i = 0
                # Insert gaps (%) = 0
                record.features.append(0) # numeric
            # Insert the other three sorting criteria: isoform type (L, B or normal), length, sp vs tr
            record.features.append(i) # numeric
            record.features.append(len(record)) # numeric
            record.features.append(record.id[:2]) # string
        #print('\n\n')
        #for record in single_species_records:
        #    print(record)
        #    print(record.features)
        #print('\n\n')
        # Sort according to all three features. Smallest comes first (ok for gaps, needs reversal for isoform and length)
        sorted_single_species_records = sorted(single_species_records, key=lambda e: (e.features[0], -e.features[1], -e.features[2], e.features[3]))
        #print(sorted_single_species_records)

        # Do we really need to evaluate "sp" (Swiss Prot) over "tr" (Translated)? It's not done in this algorhythm
        if organism_file[6:-6] not in data[3]:
            SeqIO.write(sorted_single_species_records[0], organism_file, 'fasta')

    # Sort the list according to scientific organism name
    new_records.sort(key=return_organism)
    # Write the sorted list to a new fasta file
    SeqIO.write(new_records, outfile_org, 'fasta')
    # Sort the list according to length
    new_records.sort(key=return_length)
    # Write the sorted list to a new fasta file
    SeqIO.write(new_records, outfile_len, 'fasta')

    # Assemble final fasta file for alignment and treebuilding
    all_seqs = []
    i = 0
    for file in os.listdir(output_directory):
        if not file.startswith('1_') and not file.startswith(gene):
            with open(gene+'/'+file, 'r') as readfile:
                all_seqs.append(readfile.readlines())
            i += 1
    with open(final_fasta_file, 'w') as writefile:
        for line in all_seqs:
            writefile.writelines(line)
    print_and_log_message('Remaining {0} sequences after keeping only the longest sequence if multiple sequences are available for a species: {1}:'.format(gene, i))
    for record in new_records:
        print_and_log_message('{0}:'.format(record.id))

def make_alignment(data):
    # t_coffee is too slow for large datasets! Use muscle instead.
    #execute_subprocess(
    #    "Generating multiple sequence alignment with the following command:",
    #    "t_coffee " + final_fasta_file + " -outfile " + final_aligned_fasta_file + " -output=fasta_aln -mode mcoffee")
    execute_subprocess(
        "Generating final multiple sequence alignment with the following command:",
        "muscle -in " + final_fasta_file + " -out " + final_aligned_fasta_file + " -maxiters 16")
    
def truncate_signal_peptide(data):
    # data[0] = reference sequence
    # data[4] = signal peptide length
    # data[5] = total protein length
    execute_subprocess(
        "Trimming multiple sequence alignment with the following command:",
        "t_coffee -other_pg seq_reformat -in " + final_aligned_fasta_file + " -action +extract_block \"" + data[0] + "\" " + str(data[4]) + " " + str(data[5]-1) + " > " + final_aligned_fasta_file_trimmed)

def make_tree(data):
    # File names
    final_aligned_fasta_file_trimmed = '{0}/{0}_aligned.fasta'.format(gene)
    final_aligned_fasta_file_coded = '{0}/{0}_aligned_coded.fasta'.format(gene)
    final_aligned_phylip_file_coded = '{0}/{0}_aligned_coded.phylip'.format(gene)
    final_aligned_phylip_file_decoded = '{0}/{0}_aligned_decoded_tree.phylip'.format(gene)

    execute_subprocess(
        "Converting fasta descriptions part 1 (creating code list) with t_coffee using the following command:",
        "t_coffee -other_pg seq_reformat -in " + final_aligned_fasta_file_trimmed + " -output code_name > code_names.list")

    execute_subprocess(
        "Converting fasta descriptions part 2 (replacing fasta descriptions with codes) with t_cofeee using the following command:",
        "t_coffee -other_pg seq_reformat -code code_names.list -in " + final_aligned_fasta_file_trimmed + " > " + final_aligned_fasta_file_coded)

    execute_subprocess(
        "Convert into phylip using the following command:",
        "t_coffee -other_pg seq_reformat -in " + final_aligned_fasta_file_coded + " -output phylip_aln > " + final_aligned_phylip_file_coded)

    # Detect whether parallel bootstrapping should be performed
    # -b -1: aLRT (fastest)
    # -b 1000: Bootstrap with 1000 replicates (slowest)
    mpirun_path = which('mpirun')
    phymlmpi_path = which('phyml-mpi')
    if mpirun_path != '' and phymlmpi_path != '':
        phylo_command = "mpirun -n 4 phyml-mpi -i " + final_aligned_phylip_file_coded +  " -d aa -b -1"
    else:
        phylo_command = "phyml -i " + final_aligned_phylip_file_coded +  " -d aa -b -1"

    # The gene tree building is actually never used since the species tree is used for the tree drawing.
    # We anyway calculate it to be able to compare gene and species trees.
    execute_subprocess(
        "Make tree with the following command:",
        phylo_command)

    # phyml adds or doesn't add the .txt extension to the output file (depending on the version) and we need to check for this!
    phyml_output_file = final_aligned_phylip_file_coded + "_phyml_tree"
    if os.path.isfile(phyml_output_file):
        os.rename(phyml_output_file, phyml_output_file + ".txt")
    execute_subprocess(
        "Decoding tree file file into human-readable format using the following command:",
        "t_coffee -other_pg seq_reformat -decode code_names.list -in " + phyml_output_file + ".txt > " + final_aligned_phylip_file_decoded)

def run():
    global logfile, gene, outfile_org, outfile_len, output_directory, final_fasta_file, final_aligned_fasta_file, final_aligned_fasta_file_trimmed, final_fasta_file_without_identicals
    
    # Gene dictionary:
    # referenc_sequence, min_length, max_length, blacklist, signal peptide (i.e. first aa of th emature protein), length_of_reference_sequence
    #
    # VEGF-C 33??? should be 32! check
    # VEGFB is the 186 isoform only atm!
    #
    genes = {   'VEGFA121': ['sp|P15692-9|VEGFA_HUMAN', 106, 166, [], 26, 147],
                # COMMENTS
                # tr|A0A2I3GEH5|A0A2I3GEH5_NOMLE = Northern white-cheeked gibbon, most likely bogus splicing (or a missing stop codon)
                # tr|A0A2K6MHY3|A0A2K6MHY3_RHIB = Black snub-nosed monkey, most likely bogus splicing (or a missing stop codon)
                'VEGFA165': ['sp|P15692-4|VEGFA_HUMAN', 150, 210, ['tr|A0A2I3GEH5|A0A2I3GEH5_NOMLE', 'tr|A0A2K6MHY3|A0A2K6MHY3_RHIBE'], 26, 191],
                # COMMENTS
                # tr|A0A6P5JDV6|A0A6P5JDV6_PHACI = Koala, most likely bogus splicing
                # tr|A0A6P5DSA3|A0A6P5DSA3_BOSIN = Bos indicus, most likely bogus since too many differences to the very closelu related Bos taurus
                #
                'VEGFA189': ['sp|P15692-2|VEGFA_HUMAN', 174, 234, ['tr|A0A6P5JDV6|A0A6P5JDV6_PHACI', 'tr|A0A6P5DSA3|A0A6P5DSA3_BOSIN'], 26, 215],
                'VEGFA206': ['sp|P15692|VEGFA_HUMAN', 191, 251, ['tr|Q96FD9|Q96FD9_HUMAN', 'tr|A0A674HBF7|A0A674HBF7_TAEGU'], 26, 232],
                #'PGF': ['sp|P49763-3|PLGF_HUMAN', 155, 183, ['tr|A0A667G928|A0A667G928_LYNCA', 'tr|A0A6J1ZWP5|A0A6J1ZWP5_ACIJB', 'tr|A0A6P4U8C1|A0A6P4U8C1_PANPR', 'tr|A0A673TNV8|A0A673TNV8_SURSU'], 18, 170],
                #'VEGFB': ['sp|P49765|VEGFB_HUMAN', 199, 220, ['tr|A0A6J0UH77|A0A6J0UH77_9SAUR', 'tr|A0A6J3BK97|A0A6J3BK97_VICPA', 'tr|A0A6J1YXL3|A0A6J1YXL3_ACIJB', 'tr|A0A6J3RQ20|A0A6J3RQ20_TURTR', 'tr|A0A2Y9PWA5|A0A2Y9PWA5_DELLE'], 21, 207],
                # COMMENTS
                # tr|A0A7L1ZLL3|A0A7L1ZLL3_LEILU = Red-billed leiothrix, this is only a fragment
                # tr|H0Z7S4|H0Z7S4_TAEGU = Zebra finch, massively divergent, probably wrong splice predction
                # tr|A0A7K6NMR3|A0A7K6NMR3_PEDTO = Plains-wanderer, internal deletions in VHD
                # tr|A0A7L3TLD5|A0A7L3TLD5_URIAL = Common murre, this is only a fragment
                # tr|A0A7L0WCJ7|A0A7L0WCJ7_ALELA = Australian brushturkey, this is only a fragment
                # tr|A0A384ANZ6|A0A384ANZ6_BALAS = minke whale, N-terminal massive disagreement
                # tr|A0A452VFK6|A0A452VFK6_URSMA = polar bear,N-terminal massive disagreements
                # tr|A0A1U7U4V5|A0A1U7U4V5_CARSF = Philippine tarsier, N-terminal massive disagreements 
                #
                'VEGFC': ['sp|P49767|VEGFC_HUMAN', 375, 425, ['tr|A0A7L1ZLL3|A0A7L1ZLL3_LEILU',
                    'tr|H0Z7S4|H0Z7S4_TAEGU',
                    'tr|A0A7K6NMR3|A0A7K6NMR3_PEDTO',
                    'tr|A0A7L3TLD5|A0A7L3TLD5_URIAL',
                    'tr|A0A7L0WCJ7|A0A7L0WCJ7_ALELA',
                    'tr|A0A384ANZ6|A0A384ANZ6_BALAS',
                    'tr|A0A452VFK6|A0A452VFK6_URSMA',
                    'tr|A0A1U7U4V5|A0A1U7U4V5_CARSF'], 31, 419],
                # COMMENTS
                #
                #
                # tr|A0A455BX13|A0A455BX13_PHYMC = sperm whale, fragment
                # tr|A0A7K6ZZH8|A0A7K6ZZH8_9AVES = Andean tinamou, massive frameshift or similar within VHD
                # tr|A0A7K9J7Q6|A0A7K9J7Q6_9CORV = Velvet flycatcher, N-terminal insertion
                'VEGFD': ['sp|O43915|VEGFD_HUMAN', 270, 385, ['tr|G1MFE0|G1MFE0_AILME',
                    'tr|A0A7E6D0F3|A0A7E6D0F3_9CHIR',
                    'tr|A0A7E6D135|A0A7E6D135_9CHIR',
                    'tr|A0A455BX13|A0A455BX13_PHYMC',
                    'tr|A0A7K6ZZH8|A0A7K6ZZH8_9AVES',
                    'tr|A0A7K9J7Q6|A0A7K9J7Q6_9CORV'], 21, 354]
            }
    
    now = datetime.now()
    logfile = 'analysis/analysis_{0}.log'.format(now.strftime("%Y-%m-%d_%H:%M:%S"))

    for gene, data in genes.items():
        if data != []:
            # File names
            if gene[0:5] == 'VEGFA':
                #
                # THE STUFF BELOW IS NOT APPLICABLE ANYMORE
                # Commented-out is the first run to get alignments,
                # second, active command is executed on the intersection of the
                # results from the first run (143 species, for which all A189, C and D
                # sequences are available)
                #
                infile = 'uniprot-gene_{0}_intersection.fasta'.format(gene[0:5].lower())
                #infile = '{0}_final_aligned_intersection.fasta'.format(gene[0:5].lower())
            else:
                infile = 'uniprot-gene_{0}_intersection.fasta'.format(gene.lower())
                #infile = '{0}_final_aligned_intersection.fasta'.format(gene.lower())
            outfile_org = '{0}/{0}_sorted_organism.fasta'.format(gene)
            outfile_len = '{0}/{0}_sorted_length.fasta'.format(gene)
            output_directory = '{0}'.format(gene)
            final_fasta_file = '{0}/{0}_final.fasta'.format(gene)
            final_fasta_file_without_identicals = '{0}/{0}_without_identicals.fasta'.format(gene)
            final_aligned_fasta_file = '{0}/{0}_final_aligned.fasta'.format(gene)
            final_aligned_fasta_file_trimmed = '{0}/{0}_final_aligned_noSP.fasta'.format(gene)

            clean_sequences(infile, gene, data)
            make_alignment(data)
            truncate_signal_peptide(data)
            #make_tree(data)

if __name__ == '__main__':
    run()
