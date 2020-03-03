#!/bin/bash
#$ -cwd
#$ -N aln
#$ -S /bin/bash
#$ -b y
#$ -l mem_requested=10G
#$ -l h_vmem=10G
#$ -l tmp_requested=5G
#$ -pe smp 8
# #$ -e /dev/null
# #$ -o /dev/null


module load marsmi/nanopore/minimap2/2.17-r943-dirty
module load evaben7/gcc/8.2.0
module load evaben7/samtools/1.9/gcc-8.2.0

REF=${1}
FASTQS=${2}
OUTPUT=${3}
QUERY=$(find $FASTQS -type f -name "*.fastq" | sed -n ${SGE_TASK_ID}p)
STEM=${QUERY##*/}


# build a paf and a sam
# add in better alignment for paf  - cigar string or whatever
minimap2 -x map-ont -t 8 -k15 $REF $QUERY > ${OUTPUT}/${STEM%*.fastq}.paf
minimap2 -ax map-ont -t 8 -k15 $REF $QUERY | samtools view -Sb - | samtools sort -o ${OUTPUT}/${STEM%*.fastq}.srt.bam -
samtools index ${OUTPUT}/${STEM%*.fastq}.srt.bam

# depricate this step
