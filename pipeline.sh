#!/bin/bash


# HIVE-seq

WORK_DIR=$1
RAW_DATA=$2
SAMPLE_NAME=$3
SCRIPTS=${WORK_DIR}/scripts
BAMS=${WORK_DIR}/bams
CALLS=${WORK_DIR}/calls
# TODO: make REF an arg
REF=/directflow/KCCGGenometechTemp/projects/jamfer/HIVE/ref/HXB2_trimmed.fasta
GUPPY_OUT=${WORK_DIR}/guppy_output
# TODO: Make GUP_MODEL an arg with a default
GUP_MODEL=dna_r9.4.1_450bps_hac.cfg

# TODO: extend to include all new args
if [ -z ${WORK_DIR} ]; then echo "WORK_DIR required, not present"; exit 1; fi
if [ -z ${RAW_DATA} ]; then echo "RAW_DATA required, not present"; exit 1; fi
if [ -z ${SAMPLE_NAME} ]; then echo "SAMPLE_NAME required, not present"; exit 1; fi

# TODO: remove SGE specific runs and make regular bash calls


echo -e "[SGE - $(date +"%T")]\t1. Launching basecalling"
NUM_JOBS=$(find ${RAW_DATA} -type f -name "*.fast5" | wc -l)
echo -e "[SGE - $(date +"%T")]\tNumber of jobs to launch: $NUM_JOBS"
# TODO: remove qsub
# TODO: fix bcaller.sge to be less shit
qsub -terse -sync y ${SCRIPTS}/bcaller.sge ${RAW_DATA} ${GUP_MODEL}
echo -e "[SGE - $(date +"%T")]\tSkipping basecalling for testing"

tail -1 ${GUPPY_OUT}/time.txt

echo -e "[SGE - $(date +"%T")]\t2. Merging basecalled data"
# merging fastq files
if [ ! -d ${WORK_DIR}/fastqs ]; then mkdir -p ${WORK_DIR}/fastqs; fi
FASTQS=${WORK_DIR}/fastqs
for bc in ${GUPPY_OUT}/barcode*; do $(for fq in ${bc}/*.fastq; do cat ${fq}; done >> ${FASTQS}/${bc##*/}.fastq); done
echo -e "[SGE - $(date +"%T")]\tSkipping merging for testing"
# merge seq_sum files
head -1 ${GUPPY_OUT}/sequencing_summary.txt > ${WORK_DIR}/sequencing_summary.txt
for file in ${GUPPY_OUT}/*/sequencing_summary.txt; do tail -n +2 $file; done >> ${WORK_DIR}/sequencing_summary.txt

echo -e "[SGE - $(date +"%T")]\tAligning data"
NUM_JOBS=$(find ${FASTQS} -type f -name "*.fastq" | wc -l)
echo -e "[SGE - $(date +"%T")]\tNumber of jobs to launch: $NUM_JOBS"
if [ ! -d ${BAMS} ]; then mkdir -p ${BAMS}; fi
# array job for each barcode, do alignment
# TODO: remove qsub
# TODO: fix align script
qsub -terse -sync y -t 1-${NUM_JOBS} ${SCRIPTS}/align.sh ${REF} ${FASTQS} ${BAMS}

# echo -e "[SGE - $(date +"%T")]\tSkipping aligning for testing"

echo -e "[SGE - $(date +"%T")]\tHIVE time!"
NUM_JOBS=$(find ${BAMS} -type f -name "*.bam" | wc -l)
echo -e "[SGE - $(date +"%T")]\tNumber of jobs to launch: $NUM_JOBS"
if [ ! -d ${CALLS} ]; then mkdir -p ${CALLS}; fi
# TODO: remove qsub
# TODO: Fix hive.sh
qsub -terse -sync y -t 1-${NUM_JOBS} ${SCRIPTS}/hive.sh ${BAMS} ${WORK_DIR} ${CALLS} ${FASTQS}

#TODO: change log message
echo -e "[SGE - $(date +"%T")]\tReady for next step"
