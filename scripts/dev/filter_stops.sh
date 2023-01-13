#!/bin/bash

# module load marsmi/nanopore/minimap2/2.17-r943-dirty
# module load evaben7/gcc/8.2.0
# module load evaben7/samtools/1.9/gcc-8.2.0
# module load evaben7/bcftools/1.9/gcc-8.2.0
# module load evaben7/htslib/1.9/gcc-8.2.0
#
# source /home/jamfer/work/venv363/bin/activate

# WKDIR
# ${REF}
# ${STEM}_fwdRev.Q10.8500.srt.bam
# ${STEM}_variants.csq.vcf
# ${STEM}_fwdRev.Q10.8500.fastq
# ${STEM}_with_stop.fastq
# ${STEM}_no_stop.fastq

WKDIR=$1
ROUND=$2
STEM=$3
REF=$4
BAM=$5
VAR=$6
IN_FASTQ=$7
FILTER_READS=${STEM}_filter_reads_${ROUND}.txt
WITH_STOP_FASTQ=${STEM}_with_stop_${ROUND}.fastq
NO_STOP_FASTQ=${STEM}_no_stop_${ROUND}.fastq
REGION_BED=${WKDIR}/scripts/region.bed
VARSCAN=${WKDIR}/scripts/VarScan.v2.4.3.jar


grep stop_gain ${VAR} | awk '{print $1"_"$2"_"$4"_"$5}' > ${STEM}_STOP_GAINS_THAT_WE_ARE_LOOKING_FOR_${ROUND}.txt


# for line in $(samtools view ${BAM})
samtools view ${BAM}| while read line;
do
    READ_ID="$(echo "${line}" | cut -f1)"
    samtools view -H ${BAM} > ${STEM}_tmp.sam
    echo "${line}" >> ${STEM}_tmp.sam
    samtools view -S -b ${STEM}_tmp.sam > ${STEM}_tmp_${READ_ID}.bam

    samtools index ${STEM}_tmp_${READ_ID}.bam
    samtools mpileup --min-BQ 0 -f ${REF} --max-depth 10 --positions ${REGION_BED} ${STEM}_tmp_${READ_ID}.bam -o ${STEM}_single_read.mpileup

    java -jar ${VARSCAN} mpileup2snp ${STEM}_single_read.mpileup --min-avg-qual 0 --min-coverage 0 --min-reads2 0 --min-var-freq 0.01 --p-value 0.9999 --variants --output-vcf --strand-filter 0 > ${STEM}_single_read.snvs.vcf
    java -jar ${VARSCAN} mpileup2indel ${STEM}_single_read.mpileup --min-avg-qual 0 --min-coverage 0 --min-reads2 0 --min-var-freq 0.01 --p-value 0.9999 --variants --output-vcf --strand-filter 0 > ${STEM}_single_read.indels.vcf

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${STEM}_single_read.snvs.vcf | sed 's/%//g' > ${STEM}_snv_tmp.tsv
    cat ${STEM}_snv_tmp.tsv | awk '{print $1"_"$2"_"$3"_"$4"\tSNV"}' > ${STEM}_snv_tmp2.tsv

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${STEM}_single_read.indels.vcf | sed 's/%//g' > ${STEM}_indel_tmp.tsv
    cat ${STEM}_indel_tmp.tsv | awk '{print $1"_"$2"_"$3"_"$4"\tINDEL"}' >  ${STEM}_indel_tmp2.tsv

    cat ${STEM}_snv_tmp2.tsv ${STEM}_indel_tmp2.tsv > ${STEM}_single_read.varscan_calls_${READ_ID}_${ROUND}.tsv

    python3 ${WKDIR}/scripts/test_vars.py -s ${STEM}_STOP_GAINS_THAT_WE_ARE_LOOKING_FOR_${ROUND}.txt -v ${STEM}_single_read.varscan_calls_${READ_ID}_${ROUND}.tsv -r ${READ_ID} -f ${FILTER_READS}

done

python3 ${WKDIR}/scripts/filter_stops.py -f ${IN_FASTQ} -t ${FILTER_READS} -w ${WITH_STOP_FASTQ} -n ${NO_STOP_FASTQ}

#
# samtools mpileup \
# --min-BQ 0 \
# -f HXB2_trimmed.fasta \
# --max-depth 10 \
# --positions region.bed \
# single_read.bam -o single_read.mpileup
#
# java -jar VarScan.v2.4.3.jar mpileup2snp \
# single_read.mpileup \
# --min-avg-qual 0 \
# --min-coverage 0 \
# --min-reads2 0 \
# --min-var-freq 0.01 \
# --p-value 0.9999 \
# --variants \
# --output-vcf \
# --strand-filter 0 \
# > single_read.snvs.vcf
#
# java -jar VarScan.v2.4.3.jar mpileup2indel \
# single_read.mpileup \
# --min-avg-qual 0 \
# --min-coverage 0 \
# --min-reads2 0 \
# --min-var-freq 0.01 \
# --p-value 0.9999 \
# --variants \
# --output-vcf \
# --strand-filter 0 \
# > single_read.indels.vcf
#
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' single_read.snvs.vcf | \
# sed 's/%//g' > snv_tmp.tsv
# cat snv_tmp.tsv | awk '{print $1"_"$2"_"$3"_"$4"\tSNV"}' > snv_tmp2.tsv
#
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' single_read.indels.vcf | \
# sed 's/%//g' > indel_tmp.tsv
# cat indel_tmp.tsv | awk '{print $1"_"$2"_"$3"_"$4"\tINDEL"}' >  indel_tmp2.tsv
#
# cat snv_tmp2.tsv indel_tmp2.tsv > single_read.varscan_calls.tsv
#
#
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' stop_variants.vcf | awk '{print $1"_"$2"_"$3"_"$4}' > STOP_GAINS_THAT_WE_ARE_LOOKING_FOR.txt
#
# fetchSubset.pl single_read.varscan_calls.tsv STOP_GAINS_THAT_WE_ARE_LOOKING_FOR.txt 1 1 > STOPS_THAT_ARE_PRESENT_IN_THE_READ.txt
# removeSubset.pl STOP_GAINS_THAT_WE_ARE_LOOKING_FOR.txt single_read.varscan_calls.tsv 1 1 > STOPS_THAT_ARE_NOT_PRESENT_IN_THE_READ.txt

#bcftools csq -l -pa -f HXB2_trimmed.fasta -g HXB2_AA.gff single_read.snvs.vcf -Ov -o consequence_snvs.vcf
#bcftools csq -l -pa -f HXB2_trimmed.fasta -g HXB2_AA.gff single_read.indels.vcf -Ov -o consequence_indels.vcf
