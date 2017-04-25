#! /usr/bin/env bash

#BSUB -J hisat2_example
#BSUB -o logs/hisat2_align_%J.out
#BSUB -e logs/hisat2_align_%J.err
#BSUB -R "select[mem>10] rusage[mem=10]"
#BSUB -n 12 

# change default bash paramers
# -e exit if anything fails
# -u exit if any variables are undefined
# -o exit if anything fails in a pipe
# -x print commands as they are executed

set -u -x -e -o pipefail

# define variables
dir="$HOME/dev/optional-class/rna-seq"
fastqs="$dir/raw_data/"
hisat2_idx="$dir/dbases/mm10/genome"
annotations="$dir/dbases/gencode.vM13.basic.annotation.gtf"

# download the hisat2 indexes  
# from http://ccb.jhu.edu/software/hisat2/index.shtml
# otherwise need >200Gb of RAM to generated index with genome +
# transcriptome

# wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/mm10.tar.gz


########################################
#check fastq qualities
#for fastq in "$fastqs"*.fastq.gz
#do 
#    echo "getting sequence quality info for "$fastq
#    fastqc -t 12  $fastq
#done


########################################
# align reads with hisat2 
# make directory for output (-p will not throw an error if directory
# already exists

#mkdir -p hisat2
#
#for fastq in "$fastqs"*.fastq.gz
#do 
#    echo "aligning "$fastq "with hisat2"
#    
#    # strip directory information from $fastq variable
#    outname=$(basename $fastq)
#    # strip .fastq.gz from $outname variable
#    outname=${outname/.fastq.gz/}
#
#    hisat2 -x $hisat2_idx \
#        -U $fastq \
#        -p 12 \
#        | samtools view -bS - \
#        | samtools sort - \
#        > "hisat2/"$outname".bam" 
#
#    samtools index "hisat2/"$outname".bam"
#done

########################################
# count reads overlapping exons of each gene
# use featureCounts from the subread package

mkdir -p feature_counts

featureCounts \
    -T 12 \
    -s 2 \
    -t exon \
    -g gene_id \
    -a $annotations \
    -o feature_counts/counts.txt \
    hisat2/*.bam

