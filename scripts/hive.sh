#!/bin/bash


# module load centos6.10/marsmi/nanopore/minimap2/2.17-r943-dirty
# module load centos6.10/evaben7/gcc/8.2.0
# module load centos6.10/evaben7/samtools/1.9/gcc-8.2.0
# module load centos6.10/evaben7/bcftools/1.9/gcc-8.2.0
# module load centos6.10/evaben7/htslib/1.9/gcc-8.2.0
# module load centos6.10/shacar/java/jdk-11.0.2
# module load centos6.10/jamfer/jvarkit/2b40bcc
# module load centos6.10/briglo/picard/2.9.4


# source /home/jamfer/work/venv363/bin/activate

FASTQ=$1
STEM=${FASTQ##*/}
WORK_DIR=$2
# SCRIPTS=${WORK_DIR}/scripts
BAMS=${WORK_DIR}/bams
CALLS=${WORK_DIR}/calls
REF=$3
GFF=$4
PICARD=$5

echo -e "hive.sh script is ${0}"
SCRIPT=$0
SCRIPTS=${SCRIPT%/*}
echo -e "the script path is ${SCRIPTS}"

# REF=/directflow/KCCGGenometechTemp/projects/jamfer/HIVE/ref/HXB2_trimmed.fasta
# GFF=/directflow/KCCGGenometechTemp/projects/jamfer/HIVE/ref/HXB2_AA_17Jul20.gff
echo -e "[SGE - $(date +"%T")]\tHive-seq starting"
# make bams folder
if [ ! -d ${BAMS} ]; then mkdir -p ${BAMS}; fi

echo -e "[SGE - $(date +"%T")]\tBuilding .dict file"
# java -jar picard.jar CreateSequenceDictionary R=${REF}
java -jar ${PICARD} CreateSequenceDictionary R=${REF}

echo -e "[SGE - $(date +"%T")]\tBAMS: ${BAMS}"
echo -e "[SGE - $(date +"%T")]\tWORK_DIR: ${WORK_DIR}"

minimap2 -cx map-ont --secondary=no -t 8 -k15 $REF $FASTQ > ${BAMS}/${STEM%*.fastq}.paf
minimap2 -ax map-ont --secondary=no -t 8 -k15 $REF $FASTQ | samtools view -Sb - | samtools sort -o ${BAMS}/${STEM%*.fastq}.srt.bam -
samtools index ${BAMS}/${STEM%*.fastq}.srt.bam

BAM=$(find ${BAMS} -type f -name "*.srt.bam")
echo -e "[SGE - $(date +"%T")]\tBAM: ${BAM}"
BASE=${BAM##*/}
echo -e "[SGE - $(date +"%T")]\tBASE: ${BASE}"
echo -e "[SGE - $(date +"%T")]\tFASTQ: ${FASTQS}/${BASE%.srt*}.fastq"
PAF=${BAMS}/${BASE%.srt*}.paf
echo -e "[SGE - $(date +"%T")]\tPAF: ${PAF}"
FOLD=${CALLS}/${BASE%.srt*}
echo -e "[SGE - $(date +"%T")]\tFOLD: ${FOLD}"
if [ ! -d ${FOLD} ]; then mkdir -p ${FOLD}; fi
STEM=${FOLD}/${BASE%.srt*}
echo -e "[SGE - $(date +"%T")]\tSTEM: ${STEM}"

echo -e "[SGE - $(date +"%T")]\tROUND ONE - Finding stops in alignment"


# add in first alignment step too

# Get the fwd and rev reads that span the full length
samtools view ${BAM} Human:60-80 | cut -f1 | sort -u > ${STEM}_FwdReads.txt
samtools view ${BAM} Human:8991-9001 | cut -f1 | sort -u > ${STEM}_RevReads.txt

# merge them with comm
comm -12 ${STEM}_RevReads.txt ${STEM}_FwdReads.txt > ${STEM}_FwdandRev.txt
# uniq them
cat ${STEM}_FwdandRev.txt | sort | uniq > ${STEM}_FwdandRev.uniq.txt

# pull reads from bam and filter to reads found above
samtools view ${BAM} | python ${SCRIPTS}/bam_get_reads.py -r ${STEM}_FwdandRev.uniq.txt > ${STEM}_subset.sam
# get header
samtools view -H ${BAM} > ${STEM}_Header.txt
# paste them together (this could probably be done in the python file by just printing lines starting with @ rather than skipping them)
cat ${STEM}_Header.txt ${STEM}_subset.sam > ${STEM}_subset.header.sam
# turn to bam
samtools view -S -b ${STEM}_subset.header.sam > ${STEM}_fwdRev.bam
samtools index ${STEM}_fwdRev.bam
# turn into fastq
samtools bam2fq ${STEM}_fwdRev.bam > ${STEM}_fwdRev.fastq
# filter for Q10 reads
python ${SCRIPTS}/split_qscore.py -q 10 -p ${BASE%.srt*}_fwdRev.Q10 ${STEM}_fwdRev.fastq ${FOLD}/

#Count number of reads

#TODO: also output a flat file of readID, start site, CIGAR for filtering

#count number of reads
minimap2 -cx map-ont --secondary=no -t 8 -k15 ${REF} ${STEM}_fwdRev.Q10.pass.fastq > ${STEM}_fwdRev.Q10.paf
minimap2 -ax map-ont --secondary=no -k15 -t 8 ${REF} ${STEM}_fwdRev.Q10.pass.fastq | samtools view -Sb - | samtools sort -o ${STEM}_fwdRev.Q10.srt.bam -

awk '{ if ($11>=8500) print $1}' ${STEM}_fwdRev.Q10.paf > ${STEM}_Q10.8500list.txt

# filter reads with above readlist
java -jar ${PICARD} FilterSamReads I=${STEM}_fwdRev.Q10.srt.bam O=${STEM}_fwdRev.Q10.8500.srt.pre-filter.bam READ_LIST_FILE=${STEM}_Q10.8500list.txt FILTER=includeReadList
# TODO: using paf file output, filter bam using start site and CIGAR string

# bam to sam
samtools view -h ${STEM}_fwdRev.Q10.8500.srt.pre-filter.bam > ${STEM}_fwdRev.Q10.8500.srt.pre-filter.sam
# Reads CIGAR string and removes any deletions larger than -l 250
python ${SCRIPTS}/filter_bam.py -s ${STEM}_fwdRev.Q10.8500.srt.pre-filter.sam -l 250 > ${STEM}_fwdRev.Q10.8500.srt.filtered.sam
# Sort sam to Bam
samtools sort ${STEM}_fwdRev.Q10.8500.srt.filtered.sam -T ${STEM}.tmp  > ${STEM}_fwdRev.Q10.8500.srt.filtered.bam
samtools index ${STEM}_fwdRev.Q10.8500.srt.filtered.bam
# Bam to fastq
samtools bam2fq ${STEM}_fwdRev.Q10.8500.srt.filtered.bam > ${STEM}_fwdRev.Q10.8500.fastq
# index fasta for consiquence calling
samtools faidx ${REF}
# mpileup
bcftools mpileup -Ov -f ${REF} ${STEM}_fwdRev.Q10.8500.srt.filtered.bam > ${STEM}_raw.vcf
# variant calling
bcftools call -v -Ov -m  ${STEM}_raw.vcf -o ${STEM}_raw.calls.vcf
# consiquence calling
bcftools csq -l -pa -f ${REF} -g ${GFF} ${STEM}_raw.calls.vcf -Ov -o ${STEM}_variants.csq.vcf
# get number of variants called
CALLS_TOT=$(grep ^# \
            -v ${STEM}_variants.csq.vcf -c)
echo -e "[SGE - $(date +"%T")]\tNumber of Calls: ${CALLS_TOT}"

echo -e "[SGE - $(date +"%T")]\tFiltering stops, might take a while..."
# dump stop codons detected
grep stop_gained ${STEM}_variants.csq.vcf
# for each stop, find reads that have that variant at that position and put into a list of ${STEM}_reads_to_filter_1.tsv
for POS in $(grep stop_gained ${STEM}_variants.csq.vcf | awk '{print $2}'); do echo $POS; done > ${STEM}_positions_to_check_1.tsv
# go through fastq and filter into with stop/no stop outputs (DO STOP CODON CHECK HERE!!!)
python3 ${SCRIPTS}/filter_stops2.py -f ${STEM}_fwdRev.Q10.8500.fastq -b ${STEM}_fwdRev.Q10.8500.srt.filtered.bam -t ${STEM}_positions_to_check_1.tsv -w ${STEM}_with_stop_1.fastq -n ${STEM}_no_stop_1.fastq > ${STEM}_filter_stops2_stdout_1.log


POS=''
BASE=''
# ------------------------------------------------------------------------------------

# Start again to get patient reference
echo -e "[SGE - $(date +"%T")]\tROUND TWO - Finding Patient specific reference"
# ------------------------------------------------------------------------------------

minimap2 -cx map-ont --secondary=no -t 8 -k15 ${REF} ${STEM}_no_stop_1.fastq > ${STEM}_no_stop_1.paf

PAF2=${STEM}_no_stop_1.paf

minimap2 -ax map-ont --secondary=no -k15 -t 8 ${REF} ${STEM}_no_stop_1.fastq | samtools view -Sb - | samtools sort -o ${STEM}_no_stop_1.bam -

awk '{ if ($11>=8500) print $1}' ${PAF2} > ${STEM}_no_stop_1.8500list.txt

java -jar ${PICARD} FilterSamReads I=${STEM}_no_stop_1.bam O=${STEM}_no_stop_1.srt.pre-filter.bam READ_LIST_FILE=${STEM}_no_stop_1.8500list.txt FILTER=includeReadList

samtools view -h ${STEM}_no_stop_1.srt.pre-filter.bam > ${STEM}_no_stop_1.srt.pre-filter.sam
python ${SCRIPTS}/filter_bam.py -s ${STEM}_no_stop_1.srt.pre-filter.sam -l 250 > ${STEM}_no_stop_1.srt.filtered.sam
samtools sort ${STEM}_no_stop_1.srt.filtered.sam -T ${STEM}.tmp  > ${STEM}_no_stop_1.srt.filtered.bam
samtools index ${STEM}_no_stop_1.srt.filtered.bam

# Bam to fastq
samtools bam2fq ${STEM}_no_stop_1.srt.filtered.bam > ${STEM}_no_stop_1.filtered.fastq


BAM2=${STEM}_no_stop_1.srt.filtered.bam

bcftools mpileup -Ov -f ${REF} ${BAM2} > ${STEM}_raw_2.vcf

bcftools call -v -Ov -m  ${STEM}_raw_2.vcf -o ${STEM}_raw_2.calls.vcf

bcftools csq -l -pa -f ${REF} -g ${GFF} ${STEM}_raw_2.calls.vcf -Ov -o ${STEM}_variants_2.csq.vcf

CALLS_TOT=$(grep ^# \
            -v ${STEM}_variants_2.csq.vcf -c)
echo -e "[SGE - $(date +"%T")]\tNumber of Calls: ${CALLS_TOT}"

echo -e "[SGE - $(date +"%T")]\tFiltering stops again, might take a while..."

# This checks for any stops after stops were filtered (there should be none)
grep stop_gained ${STEM}_variants_2.csq.vcf

for POS in $(grep stop_gained ${STEM}_variants_2.csq.vcf | awk '{print $2}'); do echo $POS; done > ${STEM}_positions_to_check_2.tsv

python3 ${SCRIPTS}/filter_stops2.py -f ${STEM}_no_stop_1.filtered.fastq -b ${BAM2} -t ${STEM}_positions_to_check_2.tsv -w ${STEM}_with_stop_2.fastq -n ${STEM}_no_stop_2.fastq > ${STEM}_filter_stops2_stdout_2.log

POS=''
BASE=''

STOP_READ_COUNT=$(grep ^@ ${STEM}_with_stop_2.fastq -c)

# Should be zero
echo -e "[SGE - $(date +"%T")]\tNumber of stop reads (should be zero): ${STOP_READ_COUNT}"
#

# experimental:

pysamstats -d -f ${REF} --type variation_strand ${BAM2} > ${STEM}_step_cons.txt



python ${SCRIPTS}/build_consensus.py -i ${STEM}_step_cons.txt -o ${STEM}_consensus.fa

# ------------------------------------------------------------------------------------

# Start again to get all the metrics we want
echo -e "[SGE - $(date +"%T")]\tROUND THREE - Aligning to new reference and finding stops"
# ------------------------------------------------------------------------------------

REF2=${STEM}_consensus.fa
FASTQ3=${STEM}_fwdRev.Q10.8500.fastq


minimap2 -cx map-ont --secondary=no -t 8 -k15 ${REF2} ${FASTQ3} > ${STEM}_pt_ref_mapped.paf

PAF3=${STEM}_pt_ref_mapped.paf

minimap2 -ax map-ont --secondary=no -k15 -t 8 ${REF2} ${FASTQ3} | samtools view -Sb - | samtools sort -o ${STEM}_pt_ref_mapped.bam -

awk '{ if ($11>=8500) print $1}' ${PAF3} > ${STEM}_pt_ref_mapped.8500list.txt

java -jar ${PICARD} FilterSamReads I=${STEM}_pt_ref_mapped.bam O=${STEM}_pt_ref_mapped.srt.pre-filter.bam READ_LIST_FILE=${STEM}_pt_ref_mapped.8500list.txt FILTER=includeReadList

samtools view -h ${STEM}_pt_ref_mapped.srt.pre-filter.bam > ${STEM}_pt_ref_mapped.srt.pre-filter.sam
python ${SCRIPTS}/filter_bam.py -s ${STEM}_pt_ref_mapped.srt.pre-filter.sam -l 250 > ${STEM}_pt_ref_mapped.srt.filtered.sam
samtools sort ${STEM}_pt_ref_mapped.srt.filtered.sam -T ${STEM}.tmp  > ${STEM}_pt_ref_mapped.srt.filtered.bam
samtools index ${STEM}_pt_ref_mapped.srt.filtered.bam

# Bam to fastq
samtools bam2fq ${STEM}_pt_ref_mapped.srt.filtered.bam > ${STEM}_pt_ref_mapped.filtered.fastq

BAM3=${STEM}_pt_ref_mapped.srt.filtered.bam

samtools index ${BAM3}

samtools faidx ${REF2}

bcftools mpileup -Ov -f ${REF2} ${BAM3} > ${STEM}_raw_3.vcf

bcftools call -v -Ov -m  ${STEM}_raw_3.vcf -o ${STEM}_raw_3.calls.vcf

bcftools csq -l -pa -f ${REF2} -g ${GFF} ${STEM}_raw_3.calls.vcf -Ov -o ${STEM}_variants_3.csq.vcf

CALLS_TOT=$(grep ^# \
            -v ${STEM}_variants_3.csq.vcf -c)
echo -e "[SGE - $(date +"%T")]\tNumber of Calls: ${CALLS_TOT}"

echo -e "[SGE - $(date +"%T")]\tFiltering stops again, might take a while..."

echo -e "[SGE - $(date +"%T")]\tBuilding .dict file"
java -jar ${PICARD} CreateSequenceDictionary R=${REF2}

grep stop_gained ${STEM}_variants_3.csq.vcf

for POS in $(grep stop_gained ${STEM}_variants_3.csq.vcf | awk '{print $2}'); do echo $POS; done > ${STEM}_positions_to_check_3.tsv
python3 ${SCRIPTS}/filter_stops2.py -f ${STEM}_pt_ref_mapped.filtered.fastq -b ${BAM3} -t ${STEM}_positions_to_check_3.tsv -w ${STEM}_with_stop_3.fastq -n ${STEM}_no_stop_3.fastq > ${STEM}_filter_stops2_stdout_3.log
POS=''
BASE=''


FASTQ4=${STEM}_no_stop_3.fastq

minimap2 -cx map-ont --secondary=no -t 8 -k15 ${REF2} ${FASTQ4} > ${STEM}_pt_ref_no_stop_mapped.paf

PAF4=${STEM}_pt_ref_mapped.paf

minimap2 -ax map-ont --secondary=no -k15 -t 8 ${REF2} ${FASTQ4} | samtools view -Sb - | samtools sort -o ${STEM}_pt_ref_no_stop_mapped.bam -

awk '{ if ($11>=8500) print $1}' ${PAF4} > ${STEM}_pt_ref_no_stop.8500list.txt

java -jar ${PICARD} FilterSamReads I=${STEM}_pt_ref_no_stop_mapped.bam O=${STEM}_pt_ref_no_stop_mapped.srt.pre-filter.bam READ_LIST_FILE=${STEM}_pt_ref_no_stop.8500list.txt FILTER=includeReadList

samtools view -h ${STEM}_pt_ref_no_stop_mapped.srt.pre-filter.bam > ${STEM}_pt_ref_no_stop_mapped.srt.pre-filter.sam
python ${SCRIPTS}/filter_bam.py -s ${STEM}_pt_ref_no_stop_mapped.srt.pre-filter.sam -l 250 > ${STEM}_pt_ref_no_stop_mapped.srt.filtered.sam
samtools sort ${STEM}_pt_ref_no_stop_mapped.srt.filtered.sam -T ${STEM}.tmp  > ${STEM}_pt_ref_no_stop_mapped.srt.filtered.bam
samtools index ${STEM}_pt_ref_no_stop_mapped.srt.filtered.bam

BAM4=${STEM}_pt_ref_no_stop_mapped.srt.filtered.bam


echo -e "[SGE - $(date +"%T")]\tDONE"
