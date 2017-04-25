#! /usr/bin/env bash

#BSUB -J ChIP_coverage
#BSUB -o logs/chip_coverage_%J.out
#BSUB -e logs/chip_coverage_%J.err
#BSUB -R "select[mem>4] rusage[mem=4]"
#BSUB -n 1

set -u -e -x -o pipefail

# install programs if necessary (for macOSX)
# brew install meme bowtie2 samtools 
# pip install macs2

# load modules if on Tesla
# module load meme python

# define file locations
workdir="$HOME/dev/optional-class/chip-seq"
fasta="$workdir/dbases/GRCm38.p4.genome.fa"
genome="$workdir/dbases/chrom_sizes_mm38.txt"
bams="$workdir/bams"

mkdir -p coverage

for bam in $bams/*.bam
do
    echo "getting coverage data for "$bam
    
    # extract out only fastq filename without directory
    base_name=$(basename $bam)
    out_name=${base_name/.bam}

    # make compressed bedgraph (coverage data) 
    bedtools genomecov -ibam $bam -bg \
        | gzip -c > "coverage/"$out_name".bg.gz" 

    # find enriched regions (i.e. peaks) 
    macs2 callpeak -t $bam -n $out_name --outdir "coverage/" 

    # get intervals +/- 25bp around peak summit
    bedtools slop -i "coverage/"$out_name"_summits.bed" \
        -b 25 -g $genome > "coverage/"$out_name"summits.slop.25.bed" 

    # get DNA sequences for 25bp +/- surround peak summit
    bedtools getfasta -fi $fasta \
        -bed "coverage/"$out_name"summits.slop.25.bed" \
        -fo "coverage/"$out_name"summits.slop.25.fasta"
done
