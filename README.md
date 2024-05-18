# Recommended Environments:
Both tools require the use of **bedtools**. This can be installed using conda or using a docker container which contains this code as well as other code for CDR detection.
 * conda setup (requires downloading code)
   - install: `conda install bioconda::bedtools`
   - create environment: `conda create -n bedtools bioconda::bedtools`
 * docker setup
   - run docker: `docker run --rm -it -v .:./data jmmenend:0.3.2 bash`

## strict_cdr.sh 
Processes a set of input files and generates two output files: strict_cdrs.bed and strict_transitions.bed. The strict_cdrs.bed file contains the strict CDRs (Centromere Dip Region) identified based percentile thresholds, while the strict_transitions.bed file contains the more gradual transitions in and out of CDRs.  

  `strict_cdr.sh -i <input_file> -r <hg002_merged_H1L_file> -o <output_prefix> -p <percentage> -t <transition_percentage>`  
Options:  
  -i, --input: Required. Input bed4 containing fraction_modified information from modkit pileup. 
  -r, --hg002_merged_H1L: Required. Bed Track of regions you want to find CDRs within. CDRs are identified by regions of aggregate depleted 5mC methylation.
  -o, --output_prefix: Required. A prefix for the output files. The output files will be named as <output_prefix>.strictCDR.bed and <output_prefix>.strictTransitions.bed.  
  -p, --percentage: Optional. The percentage threshold for the CDRs. Default is 10.  
  -t, --transition_percentage: Optional. The transition percentage threshold. Default is 20.  

## strict_scoring.sh
Calculate strict cdrs at each percentile. Each iteration has chance to score windows. At the end calculate CDRs based on score.
