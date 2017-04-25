#! /usr/bin/env bash

#BSUB -J find_motifs 
#BSUB -o logs/find_motifs_%J.out
#BSUB -e logs/find_motifs_%J.err
#BSUB -R "span[hosts=1] select[mem>40] rusage[mem=40]"
#BSUB -n 12

set -u -e -x -o pipefail

# install programs if necessary (for macOSX)
# brew install meme 

# load modules if on Tesla
# module load meme 

input_fasta="coverage/pol2S2_mel_chipsummits.slop.25.fasta"

mkdir -p meme


outname=$(basename $input_fasta)
outname=${outname/.fasta/}


# run MEME to find enriched motifs (will take a few minutes)
meme -dna \
    -nmotifs 1 \
    -p 12  \
    -maxsize 10000000 \
    -oc "meme/"$outname \
    $input_fasta
