#!/usr/bin/python3

# Libs...
import sys
import os

import statistics

#from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO


# Read commandline var...
search_directory=sys.argv[1]
reference_db=sys.argv[2]
file_extension=sys.argv[3]
min_taxa=int(sys.argv[4])
nthreads=int(sys.argv[5])


# Set local vars...
eval_threshold=float(1e-12)
tmp_query_file = "tmp.query.fasta"
tmp_blast_output = "tmp.blast.txt"
outfile_suffix = "cleaned.fas"

# vars for input directory processing...
Orthogroup_contigs = dict()
fasta_container = dict()

# vars for BLAST processing
contig_to_accession_map = dict()
contig_to_evalue_map = dict()



# DO STUFF...

# process the files in the input directory...
# that is:
#   build a large fasta file with ALL the contained contigs
#   build a map of orthogroups to contigs for later screening
print("processing the orthogroups from the input directory...")
count_processed=0
for file in (f for f in os.listdir(search_directory)):
    if file.endswith('.' + file_extension):
        count_processed+=1
        #print("PRUNING: ", file)
        file_to_open = os.path.join(search_directory, file)
        current_contig_ids = dict()
        for current_seq in SeqIO.parse(file_to_open, "fasta"):
            if len(current_seq)==0: continue
            current_seq.seq = current_seq.seq.ungap("-")
            fasta_container[current_seq.id]=current_seq
            current_contig_ids[current_seq.id]=1
        Orthogroup_contigs[file] = list(current_contig_ids.keys())
print("   Processed ", count_processed, " files from the input directory...\n")



# write the fasta sequences out to the temporary query file that we'll need for BLAST...
print("writing the query sequences to a file for BLAST...")
o = open(tmp_query_file, "w")
for current_seq in fasta_container.keys():
    current_seq_for_writing = fasta_container[current_seq]
    current_seq_for_writing.seq = current_seq_for_writing.seq.ungap("-")
    o.write( "{0}".format( current_seq_for_writing.format("fasta") ) )
o.close()
print("   DONE\n")



# Now, run BLAST on the query string against the database...
print("running BLAST...")
run_string = "diamond blastp -d " + reference_db + " -q " + tmp_query_file + " -o " + tmp_blast_output + " --evalue " + str(eval_threshold) + " --quiet --max-target-seqs 1 --threads " + str(nthreads)
os.system(run_string)
print("   DONE\n")



# Parse the BLAST results, build a map of the contig ID to best BLAST hit, we'll process these results next in the context of each orthogroup...
print("Parsing the BLAST results file...")
tmp_blast_output = open(tmp_blast_output, 'r')
for current_BLAST_result in tmp_blast_output.readlines():
    current_BLAST_result = current_BLAST_result.strip().split("\t")
    contig_to_accession_map[ current_BLAST_result[0] ] = current_BLAST_result[1]
    contig_to_evalue_map[ current_BLAST_result[0] ] = float( current_BLAST_result[10] )
    #print(current_BLAST_result)
    #print(current_BLAST_result[0])
    #print(current_BLAST_result[1])
    #print(current_BLAST_result[10])
    #print()
tmp_blast_output.close()
print("   DONE\n")



# Process the BLAST results for each orthogroups...
print("Processing the BLAST results for each orthogroups...")
for current_orthogroup in Orthogroup_contigs.keys():
    print("Processing BLAST results for: ", current_orthogroup, "...")
    
    Accession_set = dict()

    contig_set = Orthogroup_contigs[current_orthogroup]
    for current_contig in contig_set:
        #print("   ", current_contig)
        if current_contig not in contig_to_accession_map.keys(): continue

        species = current_contig.split("|")[0]
        #print("   ", species)

        current_accession = contig_to_accession_map[current_contig]
        #print("   ", current_accession)

        evalue_match = contig_to_evalue_map[ current_contig ] 
        #print("   ", evalue_match)

        # check to see if this accession already has matches in the orthogroup, if not make a fresh dict
        if current_accession in Accession_set.keys():
            accession_species_contig = Accession_set[current_accession]
        else: accession_species_contig = dict()

        # check to see if the species is NOT in the current set
        if species not in accession_species_contig.keys():
            accession_species_contig[species] = current_contig
        else:
            # compare the evalues for the current vs prior, keep the contig with the lowest
            prior_contig = accession_species_contig[species]
            prior_evalue = contig_to_evalue_map[ prior_contig ]
            if evalue_match < prior_evalue: accession_species_contig[species] = current_contig

        Accession_set[current_accession] = accession_species_contig

    print()

    iteration = 0
    for current_Accession in Accession_set.keys():
        print( current_Accession )
        current_species_set = Accession_set[ current_Accession ]
        print( "LEN: ", len(current_species_set.keys()))
        if len(current_species_set.keys()) >= min_taxa:
            print("   MEETS MIN TAXA!")
            iteration += 1

            # print out the cleaned file...
            basefile = current_orthogroup.split(".")[0]
            
            outfile = basefile + "_" + str(iteration) + "." + outfile_suffix
            print("writing results to: ", outfile)
            outfile = os.path.join(search_directory, outfile)
            with open(outfile, "w") as o:
                for current_species in current_species_set.keys():
                    current_contig_id = current_species_set[current_species]
                    SeqIO.write(fasta_container[current_contig_id], o, "fasta")

        else:
            print("   NO! does not meet the min number of required taxa")

    print()










sys.exit()
