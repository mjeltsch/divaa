# divaa
Re-implementation of the DIVAA software with Biopython

For our current phylogenetics manuscript, we wanted an easy and intuitive measure for amino acid diversity in multiple sequence alignments. Such has been described by Rodi et al. (Rodi, D.J.; Mandava, S.; Makowski, L. DIVAA: Analysis of Amino Acid Diversity in Multiple Aligned Protein Sequences. Bioinformatics 2004, 20, 3481–3489, doi:10.1093/bioinformatics/bth432). However, we could not find the code anywhere. Because the method is straightforward, we re-implemented it using BioPython.

The scripts 1_intersect_multiple_multiple_fastafiles.py and 2_clean_and_align.py are just here to prepare the VEGF data to run with DIVAA. Because one should not mix different VEGF-A isoforms in an alignment, we needed to come up with a method to reliably sort them. We do in in the 2_clean_and_align.py script by checking length, length minus signal peptide length, but also by comparing to a list of known VEGF-A isoforms (stored in the file VEGFA_isoform_list.py).
