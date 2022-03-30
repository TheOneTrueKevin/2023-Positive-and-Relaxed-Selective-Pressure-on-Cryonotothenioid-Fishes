#!/usr/bin/python3

###############################################################################
# LIBRARIES...
from Bio import SeqIO
import os
import sys
from pathlib import Path


# read in the command line arguments...
# 1... the directory with the pruned peptide alignements
# 2... the directory with the cds files
# 3... the extension for files we are searching for...
peptide_alignment_directory = sys.argv[1]
cds_directory = sys.argv[2]
alignment_file_extension = sys.argv[3]


# read in the CDS files from the CDS directory and store the contigs in a large dict...
cds_dict = dict()
for file in (f for f in os.listdir(cds_directory)):
    file_with_path = os.path.join(cds_directory, file)
    for current_cds in SeqIO.parse(file_with_path, "fasta"):
        cds_dict[current_cds.id]=current_cds


# go through all the peptide alignments in the directory...
for file in (f for f in os.listdir(peptide_alignment_directory)):
    path_file = Path(file)
    extension = "".join(path_file.suffixes)
    terminal_extension = Path(path_file).suffix.lstrip(".")
    #print(terminal_extension)
    last_extension_removed = str(path_file.with_suffix(""))

    #print(terminal_extension, "  ", alignment_file_extension)
    if terminal_extension != alignment_file_extension: continue
    
    cds_file = last_extension_removed + ".cds"

    # open the cds file for writing...
    cds_file = os.path.join(peptide_alignment_directory,cds_file)
    cds_file = open(cds_file,'w', newline='\n')

    # parse the peptide file, write the matching cds to the output file
    file_with_path = os.path.join(peptide_alignment_directory, file)
    for current_pep in SeqIO.parse(file_with_path, "fasta"):
        matching_cds = cds_dict[current_pep.id]
        cds_file.write("{0}".format(matching_cds.format("fasta")))
    
    # close the cds file...
    cds_file.close()


