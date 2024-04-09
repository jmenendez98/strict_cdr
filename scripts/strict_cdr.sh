#!/bin/bash

# input folder should contain ONLY methylBed files as generated using the settings specified in the README
file=$2

# hg002 merged H1L bed file path
hg002_merged_H1L=$3

# output prefix
prefix=$4

# percentage of windows to keep as CDRs
percent=$1
transition_percent=$(( percent + 10 ))
# minimum length of a CDR
min_length=4500

# make output folders
mkdir -p "windows_tmp"
mkdir -p "1kbwindow_means"
mkdir -p "temp_cdrs"
mkdir -p "strict_cdrs"
mkdir -p "transition_cdrs"

# declare variables
window_bed="windows_tmp/${prefix}.windows1000.bed"
window_bed_2="windows_tmp/${prefix}.windows1000.filter.bed"
window_mean="1kbwindow_means/${prefix}.windows1000.mean.bed"
temp_cdr_1="temp_cdrs/${prefix}.strictCDR.temp.bed"
temp_transitions="temp_cdrs/${prefix}.transition.temp.bed"
temp_transitions_1="temp_cdrs/${prefix}.transition1.temp.bed"
transitions="${prefix}.strictTransitions.bed"
strict_cdrs="${prefix}.strictCDR.bed"

echo "PREFIX: " $prefix

# 1
awk 'BEGIN { min = "unset"; max = 0 }
{
	if (min == "unset" || $2 < min) min = $2;
	if ($3 > max) max = $3;
}
END {
	# Assuming the chromosome is the same for all rows and is in the first column
	chrom = $1;
	for (i = min; i <= max; i += 1001) {
		window_start = i;
		window_end = (i + 1000 > max) ? max : i + 1000;
		print chrom "\t" window_start "\t" window_end;
	}
}' $file > $window_bed

# 2
bedtools intersect -a $window_bed \
	-b $hg002_merged_H1L | bedtools sort -i - > $window_bed_2
# 3
bedtools map -a $window_bed_2 -b $file -c 4 -o mean | \
	awk -F'\t' '$4 != "." {print}' - > $window_mean

# Reset current thresholds for each file! 
current_percent=$percent
current_transition_percent=$transition_percent

# loop for adjusting the CDR_THRESHOLD and TRANSITION_THRESHOLD 
while true; do
	# 4
	cdr_threshold=$(awk '{print $4}' $window_mean | sort -n | \
		awk -v perc=$current_percent 'BEGIN{line=-1} {all[NR]=$1} END{line=int((perc/100.0)*NR); if(line<1)line=1; print all[line]}')

	cdr_transition_threshold=$(awk '{print $4}' $window_mean | sort -n | \
		awk -v perc=$current_transition_percent 'BEGIN{line=-1} {all[NR]=$1} END{line=int((perc/100.0)*NR); if(line<1)line=1; print all[line]}')

	echo "CDR Mod Percent Threshold: ${cdr_threshold}"
	echo "Transition Mod Percent Threshold: ${cdr_transition_threshold}"
	
	# Reset min_length for each file
	current_min_length=$min_length

	# loop for adjusting the MIN_LENGTH of the CDRs
	while true; do

		# 5
		awk -v thresh=$cdr_threshold '$4 < thresh' $window_mean | bedtools merge -d 3 -i - | \
			awk -v min=$current_min_length -F'\t' 'BEGIN {FS="\t"; OFS="\t"} {if ($3-$2 > min) {print $1,$2,$3}}' - > $temp_cdr_1
		bedtools intersect -a $temp_cdr_1 -b $hg002_merged_H1L -f 1.0 | \
			awk -v min=$current_min_length -F'\t' 'BEGIN {FS="\t"; OFS="\t"} {if ($3-$2 > min) {print $1,$2,$3}}' - > $strict_cdrs

		awk -v thresh=$cdr_transition_threshold '$4 < thresh' $window_mean | \
			bedtools merge -d 3 -i - > $temp_transitions
		bedtools intersect -a $temp_transitions -b $hg002_merged_H1L -f 1.0 > $temp_transitions_1
		bedtools intersect -a $temp_transitions_1 -b $strict_cdrs -wa | \
			bedtools subtract -a - -b $strict_cdrs > $transitions
		
		if [ -s "$strict_cdrs" ]; then
			break
		elif [ $current_min_length -le 2500 ]; then
			echo "Minimum length reached its lower limit for $file. Exiting loop."
			break
		else
			echo "Strict CDRs file for $file is empty. Reducing min_length and retrying..."
			current_min_length=$(( current_min_length - 500 ))  # Reduce min_length by 500
			echo "New min_length is ${current_min_length}"
		fi
	done
	
	if [ -s "$strict_cdrs" ]; then
		if [ ! -s "$transitions" ]; then
			while true; do
				cdr_transition_threshold=$(awk '{print $4}' $window_mean | sort -n | \
					awk -v perc=$current_transition_percent 'BEGIN{line=-1} {all[NR]=$1} END{line=int((perc/100.0)*NR); if(line<1)line=1; print all[line]}')
				
				awk -v thresh=$cdr_transition_threshold '$4 < thresh' $window_mean | \
					bedtools merge -d 3 -i - > $temp_transitions
				bedtools intersect -a $temp_transitions -b $hg002_merged_H1L -f 1.0 > $temp_transitions_1
				bedtools intersect -a $temp_transitions_1 -b $strict_cdrs -wa | \
					bedtools subtract -a - -b $strict_cdrs > $transitions

				if [ -s "$transitions" ]; then
					break
				else
					echo "No Transitions, Increasing Threshold"
					current_transition_percent=$(( current_transition_percent + 5 ))
				fi
			done
		fi 
		echo "Done processing $file."
		break
	elif [ $current_percent -ge 25 ]; then
			echo "Maximum percentage reached no CDRs Detected."
			break
	else
		echo "Strict CDRs file for $file is empty. Increasing percentile thresholds..."
		current_percent=$(( current_percent + 5 ))
		current_transition_percent=$(( current_transition_percent + 5 ))
	fi

done

echo "Wrote CDRs to: ${strict_cdrs}"
echo "Wrote Transitions to: ${transitions}"

echo "Done!"
