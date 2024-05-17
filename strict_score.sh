#!/bin/bash

#set -eux -o pipefail

# Initialize variables
file=""
hg002_merged_H1L=""
prefix=""
low_percent=1
high_percent=20

# Parse command-line options
while getopts ":i:r:o:l:h:" opt; do
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
		l )
			low_percent="$OPTARG"
			;;
		h )
			high_percent="$OPTARG"
			;;
		\? )
			echo "Invalid option: $OPTARG" 1>&2
			exit 1
			;;
		: )
			echo "Invalid option: $OPTARG requires an argument" 1>&2
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
echo "prefix: $prefix"
echo "percent range: ${low_percent}-${high_percent}"

# make output folders
mkdir -p "windows"
mkdir -p "temps"

# declare output variables
window_bed="windows/${prefix}.windows1000.bed"
window_bed_2="windows/${prefix}.windows1000.filter.bed"
window_mean="windows/${prefix}.windows1000.mean.bed"
cdr_scores="${prefix}.strictScores.bed"
strict_cdrs="${prefix}.strict.bed"

# generate the 1kb windows bed from the H1L censat annotations
awk 'BEGIN { min = "unset"; max = 0 }
{
	if (min == "unset" || $2 < min) min = $2;
	if ($3 > max) max = $3;
}
END {
	# Assuming the chromosome is the same for all rows and is in the first column
	chrom = $1;
	for (i = min; i <= max; i += 500) {
		window_start = i;
		window_end = (i + 499 > max) ? max : i + 499;
		print chrom "\t" window_start "\t" window_end;
	}
}' $file > $window_bed

# generate bed
bedtools intersect -a $window_bed -b $hg002_merged_H1L | \
	sort -k 1,1 -k2,2n - | \
	bedtools map -a - -b $file -c 4 -o mean | \
	awk -F'\t' '$4 != "." {print}' - | \
	sort -k 1,1 -k2,2n - > $window_mean

# initialize cdr scoring bedgraph
awk 'OFS="\t" {print $1, $2, $3, 0}' $window_bed > $cdr_scores

# Store scoring bed file in an associative array
declare -A windows
while read -r chrom start end score; do
	windows["$chrom,$start,$end"]=$score
done < $cdr_scores

for ((percent=high_percent; percent>=low_percent; percent--)); do
	# generate the thresholds for the strict CDR/transitions from the percentages
	threshold=$(awk '{print $4}' $window_mean | sort -n | \
		awk -v perc=$percent 'BEGIN{line=-1} {all[NR]=$1} END{line=int((perc/100.0)*NR); if(line<1)line=1; print all[line]}')

	# Loop through each line of the methylation bed file
	while IFS=$'\t' read -r chrom start end methylation; do
		# Check if the average methylation is below the threshold
		if (( $(echo "$methylation < $threshold" | bc -l) )); then
			windows["$chrom,$start,$end"]=$((windows["$chrom,$start,$end"] + 1))
		fi
	done < $window_mean
done

# Write updated scoring bed file
for key in "${!windows[@]}"; do
    # Replace commas with tabs in the key
    formatted_key="${key//,/$'\t'}"
    # Echo the formatted key-value pair
    echo -e "$formatted_key\t${windows[$key]}"
done > tmpfile && mv tmpfile ${cdr_scores}

# sort output :)
sort -k 1,1 -k2,2n -o ${cdr_scores} ${cdr_scores}

# Find the Strict High Confidence CDRs (merge nearby high confidence windows with this)
awk 'BEGIN {OFS=FS="\t"} $4 >= 17 {print $1, $2, $3, $4}' ${cdr_scores} | \
	sort -k 1,1 -k2,2n - | \
	bedtools merge -d 1001 -c 4 -o 'mean' -i - | \
	awk 'BEGIN {OFS=FS="\t"} $3 - $2 >= 1750 {print $1, $2, $3, "strict_CDR", 0, ".", $2, $3, "0,0,255"}' | \
	sort -k 1,1 -k2,2n - | \
	awk '!seen[$0]++' > temp_cdrs.bed

# Find the Strict Medium Confidence CDRs (only merge adjacent windows)
awk 'BEGIN {OFS=FS="\t"} $4 > 10 {print $1, $2, $3, $4}' ${cdr_scores} | \
	sort -k 1,1 -k2,2n - | \
	bedtools merge -d 1 -c 4 -o 'mean' -i - | \
	bedtools subtract -a - -b temp_cdrs.bed | \
	awk 'BEGIN {OFS=FS="\t"} $3 - $2 >= 1750 {print $1, $2, $3, "strict_CDR", 0, ".", $2, $3, "0,0,255"}' | \
	sort -k 1,1 -k2,2n - | \
	awk '!seen[$0]++' >> temp_cdrs.bed
	
# Find the Strict Transitions
awk 'BEGIN {OFS=FS="\t"} $4 > 5 {print $1, $2, $3, $4}' ${cdr_scores} | \
	sort -k 1,1 -k2,2n - | \
	bedtools merge -d 1001 -c 4 -o 'mean' -i - | \
	bedtools intersect -a - -b temp_cdrs.bed -wa | \
	bedtools subtract -a - -b temp_cdrs.bed | \
	awk 'BEGIN {OFS=FS="\t"} {print $1, $2, $3, "strict_Transition", 0, ".", $2, $3, "173,216,230"}' | \
	sort -k 1,1 -k2,2n - | \
	awk '!seen[$0]++' > temp_transitions.bed

# generate output
cat temp_cdrs.bed temp_transitions.bed | \
	sort -k1,1 -k2,2n -o ${strict_cdrs} 