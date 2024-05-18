#!/bin/bash

set -eux -o pipefail

# Initialize variables
file=""
hg002_merged_H1L=""
prefix=""
percent="10"
transition_percent="20"

# Check if the correct number of arguments are provided
if [ "$#" -lt 1 ]; then
    echo "Error: Missing required arguments. Use -h for help."
    exit 1
fi

# Parse command-line options
while getopts ":i:r:o:p:t:" opt; do
	case ${opt} in
		i )
			file="$OPTARG"
			;;
		r )
			hg002_merged_H1L="$OPTARG"
			;;
		o )
			prefix=$(basename "$OPTARG")
			;;
		p )
			percent="$OPTARG"
			;;
		t )
			transition_percent="$OPTARG"
			;;
		h )
            echo "Usage: strict_cdr.sh -i <input_file> -r <hg002_merged_H1L_file> -o <output_prefix> -p <percentage> -t <transition_percentage>"
            echo "Options:"
            echo "  -i, --input: Required. The input file containing the CDRs and transitions."
            echo "  -r, --hg002_merged_H1L: Required. The file containing the merged H1L (Human-Mouse Homolog) regions."
            echo "  -o, --output_prefix: Required. A prefix for the output files. The output files will be named as <output_prefix>.strictCDR.bed and <output_prefix>.strictTransitions.bed."
            echo "  -p, --percentage: Optional. The percentage threshold for the CDRs. Default is 10."
            echo "  -t, --transition_percentage: Optional. The transition percentage threshold. Default is 20."
            exit 0
            ;;
		\? )
			echo "Invalid option: $OPTARG" 1>&2
			echo "Use -h for help."
			exit 1
			;;
	esac
done

# Check if required options are provided
if [[ -z "$file" || -z "$hg002_merged_H1L" || -z "$prefix" ]]; then
	echo "Missing required arguments"
	exit 1
fi

# Print parsed arguments (optional)
echo "file: $file"
echo "hg002_merged_H1L: $hg002_merged_H1L"
echo "output_prefix: $prefix"
echo "percent: $percent"
echo "transition_percent: $transition_percent"

# minimum length of a CDR
min_length=4500

# make output folders
mkdir -p "windows"
mkdir -p "temps"

# declare variables
window_bed="windows/${prefix}.windows1000.bed"
window_bed_2="windows/${prefix}.windows1000.filter.bed"
window_mean="windows/${prefix}.windows1000.mean.bed"
transitions="${prefix}.strictTransitions.bed"
strict_cdrs="${prefix}.strictCDR.bed"

# 1
awk 'BEGIN { min = "unset"; max = 0 }
{
	if (min == "unset" || $2 < min) min = $2;
	if ($3 > max) max = $3;
}
END {
	# Assuming the chromosome is the same for all rows and is in the first column
	chrom = $1;
	for (i = min; i <= max; i += 1000) {
		window_start = i;
		window_end = (i + 999 > max) ? max : i + 999;
		print chrom "\t" window_start "\t" window_end;
	}
}' $file > $window_bed

# 2 + 3
bedtools intersect -a $window_bed -b $hg002_merged_H1L | \
	sort -k 1,1 -k2,2n - | \
	bedtools map -a - -b $file -c 4 -o mean | \
	awk -F'\t' '$4 != "." {print}' - | \
	sort -k 1,1 -k2,2n - > $window_mean

# Reset current thresholds for each file! 
current_percent=$percent
current_transition_percent=$transition_percent

# 4
cdr_threshold=$(awk '{print $4}' $window_mean | sort -n | \
	awk -v perc=$current_percent 'BEGIN{line=-1} {all[NR]=$1} END{line=int((perc/100.0)*NR); if(line<1)line=1; print all[line]}')

cdr_transition_threshold=$(awk '{print $4}' $window_mean | sort -n | \
	awk -v perc=$current_transition_percent 'BEGIN{line=-1} {all[NR]=$1} END{line=int((perc/100.0)*NR); if(line<1)line=1; print all[line]}')

# Reset min_length for each file
current_min_length=$min_length

# 5
awk -v thresh=$cdr_threshold '$4 < thresh' $window_mean | \
	bedtools merge -d 3 -i - | \
	bedtools intersect -a - -b $hg002_merged_H1L -f 1.0 | \
	awk -v min=$current_min_length -F'\t' 'BEGIN {FS="\t"; OFS="\t"} {if ($3-$2 > min) {print $1,$2,$3}}' - | \
	bedtools merge -d 1500 -i - | \
	sort -k 1,1 -k2,2n -o $strict_cdrs

awk -v thresh=$cdr_transition_threshold '$4 < thresh' $window_mean | \
	bedtools merge -d 3 -i - | \
	bedtools intersect -a - -b $hg002_merged_H1L -f 1.0 | \
	bedtools intersect -a - -b $strict_cdrs -wa | \
	bedtools subtract -a - -b $strict_cdrs | \
	sort -k 1,1 -k2,2n -o $transitions

echo "Wrote CDRs to: ${strict_cdrs}"
echo "Wrote Transitions to: ${transitions}"