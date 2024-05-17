# Recommended Environments:
Both tools require the use of **bedtools**. This can be installed using conda or using a docker container which contains this code as well as other code for CDR detection.
 * conda setup (requires downloading code)
   install: `conda install bioconda::bedtools`
   create environment: `conda create -n bedtools bioconda::bedtools`
 * docker setup
   run docker: `docker run --rm -it -v .:./data jmmenend:0.3.2 bash`

## strict_cdr.sh 
Uses bedtools to find bottom percentile of windows. Requires 5 adjacent winows to be labeled as a 'strict cdr'

## strict_scoring.sh
Calculate strict cdrs at each percentile. Each iteration has chance to score windows. At the end calculate CDRs based on score.
