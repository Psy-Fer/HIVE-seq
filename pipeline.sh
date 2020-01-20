#!/bin/bash


# HIVE-seq

WORK_DIR=$1
RAW_DATA=$2
SAMPLE_NAME=$3
SCRIPTS=${WORK_DIR}/scripts
REF=${WORK_DIR}/ref/HXB2.fasta
GUPPY_OUT=${WORK_DIR}/guppy_output/
GUP_MODEL=dna_r9.4.1_450bps_hac.cfg

if [ -z ${WORK_DIR} ]; then echo "WORK_DIR required, not present"; exit 1; fi
if [ -z ${RAW_DATA} ]; then echo "RAW_DATA required, not present"; exit 1; fi
if [ -z ${SAMPLE_NAME} ]; then echo "SAMPLE_NAME required, not present"; exit 1; fi


echo -e "[SGE - $(date +"%T")]\t1. Launching basecalling"
# NUM_JOBS=$(find ${RAW_DATA} -type f -name "*.fast5" | wc -l)
# echo -e "[SGE - $(date +"%T")]\tNumber of jobs to launch: $NUM_JOBS"
qsub -terse -sync y ${SCRIPTS}/bcaller.sge ${RAW_DATA} ${GUP_MODEL}

# for f in ${GUPPY_OUT}/time.txt; do tail -1 $f; done | sort | uniq -c

echo -e "[SGE - $(date +"%T")]\t2. Merging basecalled data"
# merging fastq files
if [ ! -d ${WORK_DIR}/fastqs/ ]; then mkdir -p ${WORK_DIR}/fastqs; fi
FASTQS=${WORK_DIR}/fastqs
for bc in ${GUPPY_OUT}/barcode*; do $(for fq in ${bc}/*.fastq; do cat ${fq}; done >> ${FASTQS}/${bc##*/}.fastq); done

# merge seq_sum files
# head -1 ${GUPPY_OUT}/sequencing_summary.txt > ${WORK_DIR}/sequencing_summary.txt
# for file in ${GUPPY_OUT}/*/sequencing_summary.txt; do tail -n +2 $file; done >> ${WORK_DIR}/sequencing_summary.txt

echo -e "[SGE - $(date +"%T")]\tAligning data"
NUM_JOBS=$(find ${FASTQS} -type f -name "*.fastq" | wc -l)
echo -e "[SGE - $(date +"%T")]\tNumber of jobs to launch: $NUM_JOBS"
if [ ! -d ${WORK_DIR}/bams/ ]; then mkdir -p ${WORK_DIR}/bams; fi
# array job for each barcode, do alignment
qsub -terse -sync y -t 1-${NUM_JOBS} ${SCRIPTS}/align.sh ${REF} ${FASTQS} ${WORK_DIR}
