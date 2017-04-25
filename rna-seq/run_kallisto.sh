#! /usr/bin/env bash

# change default bash paramers
# -e exit if anything fails
# -u exit if any variables are undefined
# -o exit if anything fails in a pipe
# -x print commands as they are executed

set -u -x -e -o pipefail

# define variables
dir="$HOME/src/MOLB7621/optional-class/rna-seq"
fasta="$dir/dbases/Mus_musculus.GRCm38.rel79.cdna.all.fa.gz"
fastqs="$dir/raw_data/"

# first build the index for kallisto (takes a few minutes)

kallisto index -i $dir"/dbases/kallisto.idx" $fasta

########################################
# psuedoalign and count reads overlapping transcripts
# make directory for output (-p will not throw an error if directory
# already exists

mkdir -p kallisto

for fastq in "$fastqs"*.fastq.gz
do 
    echo "psuedoaligning "$fastq "with kallisto"
    
    # strip directory information from $fastq variable
    outname=$(basename $fastq)
    # strip .fastq.gz from $outname variable
    outname=${outname/.fastq.gz/}

    kallisto quant --single \
        -l 250 \
        -s 50 \
        -i $dir"/dbases/kallisto.idx" \
        -o "kallisto/"$outname \
        -b 5 \
         --rf-stranded \
         $fastq
done

