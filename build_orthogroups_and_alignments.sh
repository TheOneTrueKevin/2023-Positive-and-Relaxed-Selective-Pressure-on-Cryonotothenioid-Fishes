#!/bin/bash


# Read in the command line arguments needed to run ALL the programs
########################################################################

POSITIONAL_ARGS={}

while [[ $# -gt 0 ]]; do
	case $1 in
		--Orthofinder)
			ORTHOFINDER_DIR="$2" # where to find the orthofinder executable
			shift # past argument
			shift # past value
			;;
		--OriginalPeptides)
			PEPTIDE_DIR="$2" # where to find the UNPROCESSED peptide files for running orthofinder
			shift # past argument
                        shift # past value
                        ;;
		--RenamePeptides)
			RENAME_DIR="$2" # where to find processed peptide files where contigs need renaming
			shift # past argument
                        shift # past value
                        ;;
		--ProcessedPeptides)
			PROCESSED_PEPTIDE_DIR="$2" # where to find the peptide files ready for orthofinder
			shift # past argument
                        shift # past value
                        ;;
		--CDS)
			CDS_DIR="$2" # where to find the cds files that need renaming
			# these will be moved into the PROCESSED_CDS_DIR after running
			shift # past argument
                        shift # past value
                        ;;
		--ProcessedCDS)
			PROCESSED_CDS_DIR="$2" # where to find the processed cds files for building the CDS alignments
			shift # past argument
                        shift # past value
                        ;;
		--CondensingBlastDB)
			INITIAL_BLAST_DB="$2" # the BLAST db for compressing the original peptide files
			shift # past argument
                        shift # past value
                        ;;
		--PruningBlastDB)
			BLAST_DB_FOR_CONDENSING="$2" # the BLAST DB for my limited approach to paralog pruning...
			shift # past argument
                        shift # past value
                        ;;
		--NumThreads)
			NUM_THREADS="$2" # How many threads to use per analysis 
			shift # past argument
                        shift # past value
                        ;;
		--WorkingDir)
			WORKING_DIR="$2" # the directory where alignments and tree files get moved into for screening then building alignments
			shift # past argument
                        shift # past value
                        ;;
		--NumTaxa)
			NUM_TAXA="$2" # how many total species
			shift # past argument
                        shift # past value
                        ;;
		-*|--*)
			echo "Unknown Argument $1"
			exit 1
			;;
	esac
done


echo "Welcome to my Orthogroup Alignment generating scrpit, your input selections ARE..."
echo "Orthofinder is in: ${ORTHOFINDER_DIR}"
echo "Unprocessed Peptide files are in: ${PEPTIDE_DIR}"
echo "formatted but misnamed Peptide files are in: ${RENAME_DIR}"
echo "PROCESSED Peptide files are in: ${PROCESSED_PEPTIDE_DIR}"
echo "CDS files (for processing) are in: ${CDS_DIR}"
echo "CDS files (processed) are in: $(PROCESSED_CDS_DIR)"
echo "the BLAST DB for condensing the Unprocessed peptide files (usually Swissprot) is named: ${INITIAL_BLAST_DB}"
echo "the DIAMOND DB for ghetto paralog pruning is named: ${BLAST_DB_FOR_CONDENSING}"
echo "the directory for placing all the working files is: ${WORKING_DIR}"
echo "the number of threads to use is: ${NUM_THREADS}"
echo "AND... the total number of TAXA to expect for pruning processes is: ${NUM_TAXA}"
echo



########################################################################



# condense the original peptide files against swissprot (or similar) to cut redundancy in the transcriptomes...
# this also handles renaming and places the results in the PROCESSED_PEPTIDE_DIR
####################################################################
collapse_transcriptome_by_best_blast_hit.blastp.py "$PEPTIDE_DIR" "$INITIAL_BLAST_DB" fas "$PROCESSED_PEPTIDE_DIR" "$NUM_THREADS"



# rename the fasta files in the "RENAME_DIR" to get contigs appropriately named for  orthofinder / phylopypruner
######################################################################
rename_for_orthofinder.v01.py "$RENAME_DIR" fasta "$PROCESSED_PEPTIDE_DIR"



# Run Orthofinder
#######################################################################
"$ORTHOFINDER_DIR"/orthofinder -f "$PROCESSED_PEPTIDE_DIR" -M msa -t "$NUM_THREADS" -a "$NUM_THREADS" -X -z  # -z SHOULD avoid trimming the msa

# Make the working directory
mkdir -p "$WORKING_DIR"
# Establish the directory with the OrthoFinder output
for DIRNAME in "$PROCESSED_PEPTIDE_DIR"/OrthoFinder/Results_*/
do
       [ -d "${DIRNAME}" ] && dir="${DIRNAME}" && break
done

# Move Alignments to the working directory
for FILENAME in "$DIRNAME"MultipleSequenceAlignments/*.fa
do
       cp $FILENAME $WORKING_DIR
done

# Move (& rename) Trees to the working directory
for FILENAME in "$DIRNAME"/Gene_Trees/*.txt
do
FILENAME=${FILENAME##*/}
NEW_TREE_NAME="${FILENAME%_*}.tre"
       cp "$DIRNAME"/Gene_Trees/"$FILENAME" "$WORKING_DIR"/"$NEW_TREE_NAME"
done



# SANITY CHECKING
#######################################################################

# Now match up the msa and tree files. delete any alignments where there isn't a matching tree
echo "Verifying matching alignment and tree file pairs..."
for FILENAME in "$WORKING_DIR"/*.fa
do
	TREE_FILENAME="${FILENAME%.*}.tre"
	if [ ! -f "$TREE_FILENAME" ]; then
		rm "$FILENAME"
	fi
done

# ... & match in the other direction: delete any trees where there isn't a matching alignment
for FILENAME in "$WORKING_DIR"/*.tre
do
	MSA_FILENAME="${FILENAME%.*}.fa"
	if [ ! -f "$MSA_FILENAME" ]; then
		rm "$FILENAME"
	fi
done



# INITIAL PARALOG PRUNING WITH PHYLOTREEPRUNER
######################################################################
echo "PARALOG PRUNING..."

CLASS_PATH=/usr/local/bin #Location of PhyloTreePruner and Alignment Compare .class files
for FILENAME in "$WORKING_DIR"/*.fa
do
	FILENAME=${FILENAME##*/}
	ORTHOLOGY_GROUP=${FILENAME%.*}
	java -cp $CLASS_PATH PhyloTreePruner $WORKING_DIR"/"$ORTHOLOGY_GROUP".tre" $NUM_TAXA $WORKING_DIR"/"$ORTHOLOGY_GROUP".fa" 0.75 r
done

mkdir -p $WORKING_DIR"/tree_pruner"
mv "$WORKING_DIR"/*.fa_pruned.fa "$WORKING_DIR"/tree_pruner



# clean up the working directory by storing the .fa and .tre files...
######################################################################
mkdir -p $WORKING_DIR"/processed_originals"
mv "$WORKING_DIR"/*.fa "$WORKING_DIR"/*.tre $WORKING_DIR"/processed_originals"



# finish paralog pruning the tree pruner files by BLAST against Dmawsoni...
# during the pruning... export AA sequences WITHOUT gaps!!!
####################################################################
condense_alignments.v03.py "$WORKING_DIR"/tree_pruner "$BLAST_DB_FOR_CONDENSING" fa "$NUM_TAXA" "$NUM_THREADS"

#move the cleaned files into the working directory for further processing...
mv "$WORKING_DIR"/tree_pruner/*.cleaned.fas "$WORKING_DIR"
rm tmp.blast.txt tmp.query.fasta
echo ""



########################################################################################################################################



# RENAME CDS for building alignments...
echo "RENAMING CDS CONTIGS... THIS PART BARELY WORKS... IF IT FAILS HERE THINK ABOUT PREPROCESSING..."
build_matching_cds.py "$PROCESSED_PEPTIDE_DIR" fas "$CDS_DIR" fas "$PROCESSED_CDS_DIR"
echo ""



# get the CDS for the alignment...
# the CDS files will be needed for codon aligning with GUIDANCE/PRANK
#####################################################################
echo "BUILDING THE ORTHOGROUP CDS SETS..."
subset.v02.py $WORKING_DIR $PROCESSED_CDS_DIR fas
echo ""



# realign using PRANK/GUIDANCE so that we have GOOD alignmnets for judging changes in selective pressure...
####################################################################
echo "REALIGNING USING GUIDANCE/PRANK..."

mkdir -p $WORKING_DIR"/realigned"

i=0
max_jobs=15
declare -A cur_jobs=( ) # build an associative array w/ PIDs of jobs we started

for FILENAME in "$WORKING_DIR"/*.cds
do
        echo Processing "$FILENAME"...
	CLEANED_FILENAME=${FILENAME##*/}
      	echo Cleaned "$CLEANED_FILENAME"	

        BASENAME=${FILENAME##*/}
        BASENAME=${BASENAME%%.*}
        echo "$BASENAME"
	FULLPATH="$(pwd)/$WORKING_DIR"
	echo "$FULLPATH"

        if (( ${#cur_jobs[@]} >= max_jobs )); then
                wait -n # wait for at least one job to exit
                # ...and then remove any jobs that aren't running from the table
                for pid in "${!cur_jobs[@]}"; do
                        if ! kill -0 "$pid" 2>/dev/null; then
                                unset cur_jobs[$pid]
                        fi
                done
        fi

	# move the files into a temporary directory and run GUIDANCE there
	mkdir -p "$WORKING_DIR"/"$BASENAME"
	mv "$FILENAME" "$WORKING_DIR"/"$BASENAME"

	echo perl /home/iceboy/guidance.v2.02/www/Guidance/guidance.pl
	echo --seqFile ${FULLPATH}/${BASENAME}/${CLEANED_FILENAME}
	echo --msaProgram PRANK --seqType codon
	echo --outDir ${FULLPATH}/${BASENAME}
	echo --dataset ${BASENAME}_GUIDANCE
	echo --proc_num 25 --bootstraps 25
	echo ""

	run_GUIDANCE="perl /home/iceboy/guidance.v2.02/www/Guidance/guidance.pl --seqFile ${FULLPATH}/${BASENAME}/${CLEANED_FILENAME} --msaProgram PRANK --seqType codon --outDir ${FULLPATH}/${BASENAME} --dataset ${BASENAME}_GUIDANCE --proc_num 20 --bootstraps 20"
        eval "${run_GUIDANCE}" > "$FILENAME".out & cur_jobs[$!]=1

        echo ""
        echo $! : $run_GUIDANCE
        echo ""

        ((i++))

done
wait

#####################################################################################################################################3

################## CLEAN UP ############################################
for FILENAME in "$WORKING_DIR"/*.out
do
	echo Processing "$FILENAME"...
	CLEANED_FILENAME=${FILENAME##*/}
	echo Cleaned "$CLEANED_FILENAME"
	BASENAME=${CLEANED_FILENAME%%.*}

	cp "$WORKING_DIR"/"$BASENAME"/"$BASENAME"_GUIDANCE.PRANK.Without_low_SP_Col.With_Names $WORKING_DIR"/realigned"
	rm -r "$WORKING_DIR"/"$BASENAME"
done

mkdir -p "$WORKING_DIR"/Original_CDS_and_FAS
mv "$WORKING_DIR"/*.cleaned.cds "$WORKING_DIR"/Original_CDS_and_FAS
mv "$WORKING_DIR"/*.cleaned.fas "$WORKING_DIR"/Original_CDS_and_FAS

mkdir -p "$WORKING_DIR"/logs
mv "$WORKING_DIR"/*.cds.out "$WORKING_DIR"/logs






# finalize the files for HyPhy by removing any columns with gaps, verifying taxa number, and verifying min length
#####################################################################

# copy the realigned files into the working directory
for FILENAME in "$WORKING_DIR"/realigned/*.With_Names
do
       BASENAME=${FILENAME##*/}
       BASENAME=${BASENAME%%.*}
       cp $FILENAME "$WORKING_DIR"/"$BASENAME".realigned.fas
done
FinalizeForHyPhy.v03.py $WORKING_DIR realigned.fas fas 300 $NUM_TAXA 
rm "$WORKING_DIR"/*.realigned.fas
mkdir -p $WORKING_DIR"/pruned_and_aligned"
cp "$WORKING_DIR"/*.fas "$WORKING_DIR"/pruned_and_aligned












#####################################################################
echo Done
