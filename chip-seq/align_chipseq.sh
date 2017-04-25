#! /usr/bin/env bash

#BSUB -J ChIP_example
#BSUB -o logs/chip_align_%J.out
#BSUB -e logs/chip_align_%J.err
#BSUB -R "span[hosts=1] select[mem>8] rusage[mem=8]"
#BSUB -n 12

# change default bash paramers
# -e exit if anything fails
# -u exit if any variables are undefined
# -o exit if anything fails in a pipe
# -x print commands as they are executed

set -u -e -x -o pipefail

# install programs if necessary (for macOSX)
# brew install meme bowtie2 samtools fastqc 
# pip install macs2

# load modules if on Tesla
# module load samtools bowtie2 fastqc

# define file locations
workdir="$HOME/dev/optional-class/chip-seq"
fasta="$workdir/dbases/GRCm38.p4.genome.fa"
idx="$workdir/dbases/bt2.out"

# get fastq file locations
# fastq name must end in .fastq.gz
fastqs="$workdir/raw-data"

# build indicies for bowtie2 aligner
bowtie2-build -p 6 $fasta $idx 

# check fastq qualities
for fastq in $fastqs"/*.fastq.gz"
do 
    echo "checking qualites for "$fastq
    fastqc -t 12  $fastq
done

# map data and generate bam file
# -p indicates the number of CPUs to use
# increase this if you have a multicore machine
# if working on a cluster set this value to the -n paramter indicated
# in the BSUB command

mkdir -p bams/

# extract out only fastq filename without directory

for fastq in $fastqs/*.fastq.gz
do
    echo "working on "$fastq
    base_name=$(basename $fastq)
    outname=${base_name/.fastq.gz/.bam}

    bowtie2 -p 12 -U $fastq -x $idx \
      | samtools sort - \
      | samtools view -b > "bams/"$outname

    samtools index "bams/"$outname
done
