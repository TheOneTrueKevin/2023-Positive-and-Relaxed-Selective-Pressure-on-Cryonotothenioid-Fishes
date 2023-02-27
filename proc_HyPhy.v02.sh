#!/bin/bash

########################################################################
################## Command Line Arguments ##############################
########################################################################
#Set values for variables
NUM_THREADS=$1 # thread count for orthofinder, phylopypruner, and any other programs...
WORKING_DIR=$2 # the directory where alignments and tree files get moved into for screening then building alignments
OUTPUT_LABEL=$3 # the label for the output files

TREE_FILE=$4 # the TREE to use with HyPhy
HYPHY_ANALYSIS=$5 # the TYPE of HyPhy analysis
HYPHY_ADENDUM=$6 # the trailing string, any extras for the analysis 
########################################################################


########################################################################
################## RUN HYPHY ###########################################
i=0
max_jobs=$NUM_THREADS
declare -A cur_jobs=( ) # build an associative array w/ PIDs of jobs we started
for FILENAME in "$WORKING_DIR"/*.fas
do
        FILENAME=${FILENAME##*/}
        echo Processing "$FILENAME"...
        if (( ${#cur_jobs[@]} >= max_jobs )); then
                wait -n # wait for at least one job to exit
                # ...and then remove any jobs that aren't running from the table
                for pid in "${!cur_jobs[@]}"; do
                        if ! kill -0 "$pid" 2>/dev/null; then
                                unset cur_jobs[$pid]
                        fi
                done
        fi
	HyPhy_CMD="hyphy $HYPHY_ANALYSIS --alignment $WORKING_DIR/$FILENAME --tree $WORKING_DIR/$TREE_FILE $HYPHY_ADENDUM"
	eval "${HyPhy_CMD}" > "$FILENAME".out & cur_jobs[$!]=1
	echo ""
	echo $! : $HyPhy_CMD
	echo ""

	#hyphy "$HYPHY_ANALYSIS" --alignment "$WORKING_DIR"/"$FILENAME" --tree "$WORKING_DIR"/"$TREE_FILE" --test Foreground --reference Background > "$FILENAME".out & cur_jobs[$!]=1
        #hyphy "$HYPHY_ANALYSIS" --alignment "$WORKING_DIR"/"$FILENAME" --tree "$WORKING_DIR"/"$TREE_FILE" "$HYPHY_ADENDUM" > "$FILENAME".out & cur_jobs[$!]=1
        ((i++))
done
wait


########################################################################
################## CLEAN UP ############################################


# rename the jason files...
for FILENAME in "$WORKING_DIR"/*.json
do
	FILENAME=${FILENAME##*/}
	BASENAME=${FILENAME%%.*}
	#echo "$FILENAME $BASENAME"
	mv "$WORKING_DIR"/"$FILENAME" "$WORKING_DIR"/"$BASENAME"."$OUTPUT_LABEL".json
done
mkdir -p $WORKING_DIR"/json_"$OUTPUT_LABEL
mv $WORKING_DIR"/"*.json $WORKING_DIR"/json_"$OUTPUT_LABEL


# rename and move the out files...
for FILENAME in *.out
do
       FILENAME=${FILENAME##*/}
       BASENAME=${FILENAME%%.*}
       #echo "$FILENAME $BASENAME"
       mv $FILENAME $WORKING_DIR"/json_"$OUTPUT_LABEL"/"$BASENAME"."$OUTPUT_LABEL".out"
done



# ALL DONE ;-)
#####################################################################
echo Done

