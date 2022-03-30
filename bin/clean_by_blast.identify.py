#!/usr/bin/python3

# Libs...
import sys
import os
import shutil

from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import SeqIO


# Read commandline var...
search_directory=sys.argv[1]
reference_db=sys.argv[2]
file_extension=sys.argv[3]
out_dir=sys.argv[4]


# Set local vars...
mismatches = dict()
eval_threshold=float(1e-8)
tmp_blast_output = "tmp.blastx.txt"
nthreads=20
#reference_db = "Dmaw.fas"
#tst_file = "OG0000073.fas"

# Do stuff...
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

count_mismatched=0
count_processed=0
count_empty=0
for file in (f for f in os.listdir(search_directory)):
    if file.endswith('.' + file_extension):

        count_processed+=1
        orthogroup = file.split(".")[0]
        file = os.path.join(search_directory, file)

        # sanity check, is the file empty?
        is_empty=0
        for current_seq in SeqIO.parse(file, "fasta"):
            if len(current_seq)==0:
                is_empty=1
                break
        count_empty+=is_empty
        if is_empty == 1:
            count_empty+=1
            shutil.move(file,out_dir)
            continue

        blastx_cmd = NcbiblastxCommandline(query=file, db=reference_db, evalue=eval_threshold,outfmt=6, max_target_seqs=1, out=tmp_blast_output, num_threads=nthreads)
        blastx_cmd()

        current_match=""
        with open(tmp_blast_output, 'r') as f:
            for line in f:
                line=line.strip().split("\t")
                if current_match=="": current_match=line[1]
                if current_match!=line[1]:
                    mismatches[orthogroup]=1
                    count_mismatched+=1
                    shutil.move(file,out_dir)
                    break
        os.remove(tmp_blast_output)
        print("PROCESSED:  ",count_processed)
        print("MISMATCHED: ",count_mismatched)
        print("EMPTY: ",count_empty)

# 4- report the mismatched sequences
print("PROCESSED:  ",count_processed)
print("MISMATCHED: ",count_mismatched)
print("EMPTY: ",count_empty)

# 5- remove those files with mismatches


print("All Done, bye-bye ;-)")

