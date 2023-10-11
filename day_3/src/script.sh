#!/bin/bash


# we define the input parameters here using the $ and consecutive numbers
# we will call this script like: script.sh my_reference.fasta my_fastq1.fastq my_fastq2.fastq my_results
reference=$1
fastq1=$2
fastq2=$3
output=$4


# TODO: here you need to put in the script both the alignment and the variant calling commands
# you will need to use the input parameters rather that write in this file fixed file paths or names
# this will allow you to reuse this script over any sample

# script.sh reference.fasta fastq1.fastq fastq2.fastq mutations

# we call the alignment using $reference, $fastq1 and $fastq2
# TODO: COMPLETE with the alignment command

# we call the variant calling using $reference and $output
# TODO: COMPLETE with the variant calling command
