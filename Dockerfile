# Use bedtools as the base image
FROM ubuntu:24.04

MAINTAINER Julian Menendez, jmmenend@ucsc.edu

# install bedtools v2.31.1
RUN apt-get update && \
	apt-get install bedtools=2.31.1+dfsg-2

# Copy your strict CDR script and make it accessible!
WORKDIR /opt/
COPY scripts/strict_cdr.sh .
# Make the script executable
RUN chmod +x ./strict_cdr.sh

ENV PATH=/opt:$PATH

WORKDIR /data