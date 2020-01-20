#!/bin/bash
#$ -cwd
#$ -N HIVE
#$ -S /bin/bash
#$ -b y
#$ -l mem_requested=8G
#$ -l h_vmem=8G
#$ -l tmp_requested=5G
#$ -pe smp 32
# #$ -e /dev/null
# #$ -o /dev/null

module load marsmi/nanopore/minimap2/2.17-r943-dirty
module load evaben7/gcc/8.2.0
module load evaben7/samtools/1.9/gcc-8.2.0
module load evaben7/bcftools/1.9/gcc-8.2.0

source /home/jamfer/work/venv2714/bin/activate

BAMS=$1
WORK_DIR=$2
REF=/directflow/KCCGGenometechTemp/projects/jamfer/HIVE/ref/HXB2_trimmed.fasta
GFF=/directflow/KCCGGenometechTemp/projects/jamfer/HIVE/ref/HXB2_AA.gff

BAM=$(find ${BAMS} -type f -name "*.srt.bam" | sed -n ${SGE_TASK_ID}p)
BASE=${BAM##*/}
PAF=${BAMS}/${BASE%.srt*}.paf
FOLD=${WORK_DIR}/bams/${BASE%.srt*}
if [ ! -d ${FOLD} ]; then mkdir -p ${FOLD}; fi
STEM=${FOLD}/${BASE%.srt*}

samtools view ${BAM} Human:60-80 | cut -f1 | sort -u > ${STEM}_FwdReads.txt

samtools view ${BAM} Human:8991-9001 | cut -f1 | sort -u > ${STEM}_RevReads.txt

comm -12 ${STEM}_RevReads.txt ${STEM}_FwdReads.txt > ${STEM}_FwdandRev.txt

cat ${STEM}_FwdandRev.txt | sort | uniq > ${STEM}_FwdandRev.uniq.txt


samtools view ${BAM} | python /directflow/KCCGGenometechTemp/projects/jamfer/tools/bam_get_reads.py -r ${STEM}_FwdandRev.uniq.txt > ${STEM}_subset.sam

samtools view -H ${BAM} > ${STEM}_Header.txt

cat ${STEM}_Header.txt ${STEM}_subset.sam > ${STEM}_subset.header.sam

samtools view -S -b ${STEM}_subset.header.sam > ${STEM}_fwdRev.bam

samtools index ${STEM}_fwdRev.bam

samtools bam2fq ${STEM}_fwdRev.bam > ${STEM}_fwdRev.fastq

# nanofilt -q 10 barcode05.pass.HXB2.fwdRev.fastq > barcode05.pass.HXB2.fwdRev.Q10.fastq
python /directflow/KCCGGenometechTemp/projects/jamfer/scripts/read_qbin_filt.py -f ${STEM}_fwdRev.fastq -s ${WORK_DIR}/guppy_output/sequencing_summary.txt -q 10.0 > ${STEM}_fwdRev.Q10.fastq

#Count number of reads
# make sure total map over 8500 using paf
# use paf dic here, some quality control too. Q10 + map scores?
# get on the sauce and smash it with some techno
# nanofilt -l 8500 ${STEM}_fwdRev.q10.fastq > ${STEM}_fwdRev.Q10.8500.fastq
python ${WORK_DIR}/scripts/length_paf.py -p ${PAF} -f ${STEM}_fwdRev.Q10.fastq -l 8500 >  ${STEM}_fwdRev.Q10.8500.fastq

#count number of reads

minimap2 -ax map-ont -k15 -t 8 ${REF} ${STEM}_fwdRev.Q10.8500.fastq | samtools view -Sb - | samtools sort -o ${STEM}_fwdRev.Q10.8500.srt.bam -

# samtools sort ${STEM}_fwdRev.Q10.8500.bam -T ${STEM}.tmp  > ${STEM}_fwdRev.Q10.8500.srt.bam

samtools index ${STEM}_fwdRev.Q10.8500.srt.bam


samtools faidx ${REF}

bcftools mpileup -Ov -f ${REF} ${STEM}_fwdRev.Q10.8500.srt.bam > ${STEM}_raw.vcf

bcftools call -v -Ov -m  ${STEM}_raw.vcf -o ${STEM}_raw.calls.vcf

bcftools csq -pa -f ${REF} -g ${GFF} ${STEM}_raw.calls.vcf -Ov -o ${STEM}_variants.csq.vcf



# R commands to get deletion site mapping points

#
# NanoFilt -q 10
