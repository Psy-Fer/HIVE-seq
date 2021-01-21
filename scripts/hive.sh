#!/bin/bash
#$ -cwd
#$ -N HIVE
#$ -S /bin/bash
#$ -b y
#$ -l mem_requested=16G
#$ -l h_vmem=16G
#$ -l tmp_requested=5G
#$ -pe smp 8
# #$ -e /dev/null
# #$ -o /dev/null

module load centos6.10/marsmi/nanopore/minimap2/2.17-r943-dirty
module load centos6.10/evaben7/gcc/8.2.0
module load centos6.10/evaben7/samtools/1.9/gcc-8.2.0
module load centos6.10/evaben7/bcftools/1.9/gcc-8.2.0
module load centos6.10/evaben7/htslib/1.9/gcc-8.2.0
module load centos6.10/shacar/java/jdk-11.0.2
module load centos6.10/jamfer/jvarkit/2b40bcc
module load centos6.10/briglo/picard/2.9.4


# samtools view -h barcode09_fwdRev.Q10.8500.srt.bam | sam2tsv.jar -r ../../ref/HXB2_trimmed.fasta | awk '($8==5435) && ($6=="T")' > reads_with_T_5435.tsv



source /home/jamfer/work/venv363/bin/activate

BAMS=$1
WORK_DIR=$2
CALLS=$3
FASTQS=$4
REF=/directflow/KCCGGenometechTemp/projects/jamfer/HIVE/ref/HXB2_trimmed.fasta
GFF=/directflow/KCCGGenometechTemp/projects/jamfer/HIVE/ref/HXB2_AA_17Jul20.gff
echo -e "[SGE - $(date +"%T")]\tHive-seq starting"

echo -e "[SGE - $(date +"%T")]\tBuilding .dict file"
java -jar /share/ClusterShare/software/contrib/briglo/picard/build/libs/picard.jar CreateSequenceDictionary R=${REF}

echo -e "[SGE - $(date +"%T")]\tBAMS: ${BAMS}"
echo -e "[SGE - $(date +"%T")]\tWORK_DIR: ${WORK_DIR}"


BAM=$(find ${BAMS} -type f -name "*.srt.bam" | sed -n ${SGE_TASK_ID}p)
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


samtools view ${BAM} Human:60-80 | cut -f1 | sort -u > ${STEM}_FwdReads.txt

samtools view ${BAM} Human:8991-9001 | cut -f1 | sort -u > ${STEM}_RevReads.txt

comm -12 ${STEM}_RevReads.txt ${STEM}_FwdReads.txt > ${STEM}_FwdandRev.txt

cat ${STEM}_FwdandRev.txt | sort | uniq > ${STEM}_FwdandRev.uniq.txt


samtools view ${BAM} | python ${WORK_DIR}/scripts/bam_get_reads.py -r ${STEM}_FwdandRev.uniq.txt > ${STEM}_subset.sam

samtools view -H ${BAM} > ${STEM}_Header.txt

cat ${STEM}_Header.txt ${STEM}_subset.sam > ${STEM}_subset.header.sam

samtools view -S -b ${STEM}_subset.header.sam > ${STEM}_fwdRev.bam

samtools index ${STEM}_fwdRev.bam

samtools bam2fq ${STEM}_fwdRev.bam > ${STEM}_fwdRev.fastq

# nanofilt -q 10 barcode05.pass.HXB2.fwdRev.fastq > barcode05.pass.HXB2.fwdRev.Q10.fastq
python ${WORK_DIR}/scripts/qfilter.py -f ${STEM}_fwdRev.fastq -s ${WORK_DIR}/guppy_output/sequencing_summary.txt -q 10.0 > ${STEM}_fwdRev.Q10.fastq

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

bcftools csq -l -pa -f ${REF} -g ${GFF} ${STEM}_raw.calls.vcf -Ov -o ${STEM}_variants.csq.vcf

CALLS_TOT=$(grep ^# -v ${STEM}_variants.csq.vcf -c)
echo -e "[SGE - $(date +"%T")]\tNumber of Calls: ${CALLS_TOT}"

#for f in ./*/*variants.csq.vcf; do echo $f; grep stop_gained $f | cut -f 2,4,5; done

# grep stop_gained ${STEM}_variants.csq.vcf > ${STEM}_stops_gained.tsv
#
# SGCOUNT=$(wc -l ${STEM}_stops_gained.tsv)
# echo -e "[SGE - $(date +"%T")]\tstops_gained count: ${SGCOUNT}"

#bedops --element-of 1 <(bam2bed < reads.bam) variants.bed > overlapping_reads.bed
#bedops --element-of 1 <(bam2bed < reads.bam) <(vcf2bed < variants.vcf) > overlapping_reads.bed

echo -e "[SGE - $(date +"%T")]\tFiltering stops, might take a while..."
grep stop_gained ${STEM}_variants.csq.vcf

for POS in $(grep stop_gained ${STEM}_variants.csq.vcf | awk '{print $2}'); do BASE=$(grep ${POS} ${STEM}_variants.csq.vcf | awk '{print $5}'); samtools view -h ${STEM}_fwdRev.Q10.8500.srt.bam | sam2tsv.jar -r ${REF} | awk -v a=${POS} -v b=${BASE} '($8==a) && ($6==b)'; done > ${STEM}_reads_to_filter_1.tsv
python3 ${WORK_DIR}/scripts/filter_stops.py -f ${STEM}_fwdRev.Q10.8500.fastq -t ${STEM}_reads_to_filter_1.tsv -w ${STEM}_with_stop_1.fastq -n ${STEM}_no_stop_1.fastq

POS=''
BASE=''

# grep stop_gained ${STEM}_variants.csq.vcf | python3 ${WORK_DIR}/scripts/filter_stops.py ${STEM}_fwdRev.Q10.8500.srt.bam ${STEM}_fwdRev.Q10.8500.fastq  ${STEM}_with_stop.fastq > ${STEM}_no_stop.fastq
# ${WORK_DIR}/scripts/filter_stops.sh ${WORK_DIR} 1 ${STEM} ${REF} ${STEM}_fwdRev.Q10.8500.srt.bam ${STEM}_variants.csq.vcf ${STEM}_fwdRev.Q10.8500.fastq
# ${WORK_DIR}/scripts/filter_stops.sh ${WORK_DIR} 1 ${STEM} ${REF} ${STEM}_fwdRev.Q10.8500.srt.bam ${STEM}_variants.csq.vcf ${STEM}_fwdRev.Q10.8500.fastq

# ------------------------------------------------------------------------------------

# Start again to get patient reference
echo -e "[SGE - $(date +"%T")]\tROUND TWO - Finding Patient specific reference"
# ------------------------------------------------------------------------------------

minimap2 -cx map-ont -t 8 -k15 ${REF} ${STEM}_no_stop_1.fastq > ${STEM}_no_stop_1.paf
minimap2 -ax map-ont -k15 -t 8 ${REF} ${STEM}_no_stop_1.fastq | samtools view -Sb - | samtools sort -o ${STEM}_no_stop_1.srt.bam -


BAM2=${STEM}_no_stop_1.srt.bam
samtools index ${BAM2}
PAF2=${STEM}_no_stop_1.paf

# samtools view ${BAM2} Human:60-80 | cut -f1 | sort -u > ${STEM}_FwdReads_2.txt
#
# samtools view ${BAM2} Human:8991-9001 | cut -f1 | sort -u > ${STEM}_RevReads_2.txt
#
# comm -12 ${STEM}_RevReads_2.txt ${STEM}_FwdReads_2.txt > ${STEM}_FwdandRev_2.txt
#
# cat ${STEM}_FwdandRev_2.txt | sort | uniq > ${STEM}_FwdandRev.uniq_2.txt
#
#
# samtools view ${BAM2} | python ${WORK_DIR}/scripts/bam_get_reads.py -r ${STEM}_FwdandRev.uniq_2.txt > ${STEM}_subset_2.sam
#
# samtools view -H ${BAM2} > ${STEM}_Header_2.txt
#
# cat ${STEM}_Header_2.txt ${STEM}_subset_2.sam > ${STEM}_subset.header_2.sam
#
# samtools view -S -b ${STEM}_subset.header_2.sam > ${STEM}_fwdRev_2.bam
#
# samtools index ${STEM}_fwdRev_2.bam
#
# samtools bam2fq ${STEM}_fwdRev_2.bam > ${STEM}_fwdRev_2.fastq
#
# # nanofilt -q 10 barcode05.pass.HXB2.fwdRev.fastq > barcode05.pass.HXB2.fwdRev.Q10.fastq
# python ${WORK_DIR}/scripts/qfilter.py -f ${STEM}_fwdRev_2.fastq -s ${WORK_DIR}/guppy_output/sequencing_summary.txt -q 10.0 > ${STEM}_fwdRev.Q10_2.fastq
#
# #Count number of reads
# # make sure total map over 8500 using paf
# # use paf dic here, some quality control too. Q10 + map scores?
# # get on the sauce and smash it with some techno
# # nanofilt -l 8500 ${STEM}_fwdRev.q10.fastq > ${STEM}_fwdRev.Q10.8500.fastq
# python ${WORK_DIR}/scripts/length_paf.py -p ${PAF2} -f ${STEM}_fwdRev.Q10_2.fastq -l 8500 >  ${STEM}_fwdRev.Q10.8500_2.fastq
#
# #count number of reads
#
# minimap2 -ax map-ont -k15 -t 8 ${REF} ${STEM}_fwdRev.Q10.8500_2.fastq | samtools view -Sb - | samtools sort -o ${STEM}_fwdRev.Q10.8500_2.srt.bam -
#
# # samtools sort ${STEM}_fwdRev.Q10.8500.bam -T ${STEM}.tmp  > ${STEM}_fwdRev.Q10.8500.srt.bam
#
# samtools index ${STEM}_fwdRev.Q10.8500_2.srt.bam

# samtools faidx ${REF}

# bcftools mpileup -Ov -f ${REF} ${STEM}_fwdRev.Q10.8500_2.srt.bam > ${STEM}_raw_2.vcf
bcftools mpileup -Ov -f ${REF} ${BAM2} > ${STEM}_raw_2.vcf

bcftools call -v -Ov -m  ${STEM}_raw_2.vcf -o ${STEM}_raw_2.calls.vcf

bcftools csq -l -pa -f ${REF} -g ${GFF} ${STEM}_raw_2.calls.vcf -Ov -o ${STEM}_variants_2.csq.vcf

CALLS_TOT=$(grep ^# -v ${STEM}_variants_2.csq.vcf -c)
echo -e "[SGE - $(date +"%T")]\tNumber of Calls: ${CALLS_TOT}"

echo -e "[SGE - $(date +"%T")]\tFiltering stops again, might take a while..."

# This checks for any stops after stops were filtered (there should be none)

# grep stop_gained ${STEM}_variants_2.csq.vcf | python3 ${WORK_DIR}/scripts/filter_stops.py ${STEM}_fwdRev.Q10.8500_2.srt.bam ${STEM}_no_stop.fastq  ${STEM}_with_stop_2.fastq > ${STEM}_no_stop_2.fastq
# grep stop_gained ${STEM}_variants_2.csq.vcf | python3 ${WORK_DIR}/scripts/filter_stops.py ${BAM2} ${STEM}_no_stop_1.fastq  ${STEM}_with_stop_2.fastq > ${STEM}_no_stop_2.fastq
grep stop_gained ${STEM}_variants_2.csq.vcf

for POS in $(grep stop_gained ${STEM}_variants_2.csq.vcf | awk '{print $2}'); do BASE=$(grep ${POS} ${STEM}_variants_2.csq.vcf | awk '{print $5}'); samtools view -h ${BAM2} | sam2tsv.jar -r ${REF} | awk -v a=${POS} -v b=${BASE} '($8==a) && ($6==b)'; done > ${STEM}_reads_to_filter_2.tsv
python3 ${WORK_DIR}/scripts/filter_stops.py -f ${STEM}_no_stop_1.fastq -t ${STEM}_reads_to_filter_2.tsv -w ${STEM}_with_stop_2.fastq -n ${STEM}_no_stop_2.fastq

POS=''
BASE=''
# ${WORK_DIR}/scripts/filter_stops.sh ${WORK_DIR} 2 ${STEM} ${REF} ${BAM2} ${STEM}_variants_2.csq.vcf ${STEM}_no_stop_1.fastq

STOP_READ_COUNT=$(grep ^@ ${STEM}_with_stop_2.fastq -c)

# Should be zero
echo -e "[SGE - $(date +"%T")]\tNumber of stop reads (should be zero): ${STOP_READ_COUNT}"
#
# # bcftools index ${STEM}_raw_2.calls.vcf
# bgzip -c ${STEM}_raw_2.calls.vcf > ${STEM}_raw_2.calls.vcf.gz
#
# # normalize indels
# bcftools norm -f ${REF} ${STEM}_raw_2.calls.vcf.gz -Ob -o ${STEM}_raw_2.calls.norm.bcf
#
# # filter adjacent indels within 5bp
# bcftools filter --IndelGap 5 ${STEM}_raw_2.calls.norm.bcf -Ob -o ${STEM}_raw_2.calls.norm.flt-indels.bcf
#
# bcftools index ${STEM}_raw_2.calls.norm.flt-indels.bcf
#
# cat ${REF} | bcftools consensus ${STEM}_raw_2.calls.norm.flt-indels.bcf > ${STEM}_consensus.fa

# experimental:

pysamstats -d -f ${REF} --type variation_strand ${BAM2} > ${STEM}_step_cons.txt

# For each row if column 4 value is less than 10 = TRUE enter the value in column 3
# If column 34>column 40 > column 46 > column 52 = TRUE enter A
# If column 40 > column 34 & 46 & 52 = TRUE enter C
# If column 46 > column 34 & 40 & 52 = TRUE enter T
# If column 52 > column 34 & 40 & 46 = TRUE enter G

python ${WORK_DIR}/scripts/build_consensus.py -i ${STEM}_step_cons.txt -o ${STEM}_consensus.fa

# ------------------------------------------------------------------------------------

# Start again to get all the metrics we want
echo -e "[SGE - $(date +"%T")]\tROUND THREE - Aligning to new reference and finding stops"
# ------------------------------------------------------------------------------------

REF2=${STEM}_consensus.fa
FASTQ3=${STEM}_fwdRev.Q10.8500.fastq

# minimap2 -x map-ont -t 8 -k15 ${REF2} ${FASTQS}/${BASE%.srt*}.fastq > ${STEM}_no_stop_2.paf
# minimap2 -ax map-ont -k15 -t 8 ${REF2} ${FASTQS}/${BASE%.srt*}.fastq | samtools view -Sb - | samtools sort -o ${STEM}_no_stop_2.srt.bam -

minimap2 -cx map-ont -t 8 -k15 ${REF2} ${FASTQ3} > ${STEM}_pt_ref_mapped.paf
minimap2 -ax map-ont -k15 -t 8 ${REF2} ${FASTQ3} | samtools view -Sb - | samtools sort -o ${STEM}_pt_ref_mapped.srt.bam -

BAM3=${STEM}_pt_ref_mapped.srt.bam
samtools index ${BAM3}
PAF3=${STEM}_pt_ref_mapped.paf

# samtools view ${BAM2} Human:60-80 | cut -f1 | sort -u > ${STEM}_FwdReads_3.txt
#
# samtools view ${BAM2} Human:8991-9001 | cut -f1 | sort -u > ${STEM}_RevReads_3.txt
#
# comm -12 ${STEM}_RevReads_3.txt ${STEM}_FwdReads_3.txt > ${STEM}_FwdandRev_3.txt
#
# cat ${STEM}_FwdandRev_3.txt | sort | uniq > ${STEM}_FwdandRev.uniq_3.txt
#
#
# samtools view ${BAM2} | python ${WORK_DIR}/scripts/bam_get_reads.py -r ${STEM}_FwdandRev.uniq_3.txt > ${STEM}_subset_3.sam
#
# samtools view -H ${BAM2} > ${STEM}_Header_3.txt
#
# cat ${STEM}_Header_3.txt ${STEM}_subset_3.sam > ${STEM}_subset.header_3.sam
#
# samtools view -S -b ${STEM}_subset.header_3.sam > ${STEM}_fwdRev_3.bam
#
# samtools index ${STEM}_fwdRev_3.bam
#
# samtools bam2fq ${STEM}_fwdRev_3.bam > ${STEM}_fwdRev_3.fastq
#
#
# python ${WORK_DIR}/scripts/qfilter.py -f ${STEM}_fwdRev_3.fastq -s ${WORK_DIR}/guppy_output/sequencing_summary.txt -q 10.0 > ${STEM}_fwdRev.Q10_3.fastq
#
#
# python ${WORK_DIR}/scripts/length_paf.py -p ${PAF3} -f ${STEM}_fwdRev.Q10_3.fastq -l 8500 >  ${STEM}_fwdRev.Q10.8500_3.fastq
#
#
# minimap2 -ax map-ont -k15 -t 8 ${REF2} ${STEM}_fwdRev.Q10.8500_3.fastq | samtools view -Sb - | samtools sort -o ${STEM}_fwdRev.Q10.8500_3.srt.bam -
#
#
samtools index ${BAM3}

samtools faidx ${REF2}

bcftools mpileup -Ov -f ${REF2} ${BAM3} > ${STEM}_raw_3.vcf

bcftools call -v -Ov -m  ${STEM}_raw_3.vcf -o ${STEM}_raw_3.calls.vcf

bcftools csq -l -pa -f ${REF2} -g ${GFF} ${STEM}_raw_3.calls.vcf -Ov -o ${STEM}_variants_3.csq.vcf

CALLS_TOT=$(grep ^# -v ${STEM}_variants_3.csq.vcf -c)
echo -e "[SGE - $(date +"%T")]\tNumber of Calls: ${CALLS_TOT}"

echo -e "[SGE - $(date +"%T")]\tFiltering stops again, might take a while..."
# grep stop_gained ${STEM}_variants_3.csq.vcf | python3 ${WORK_DIR}/scripts/filter_stops.py ${BAM3} ${FASTQ3} ${STEM}_with_stop_3.fastq > ${STEM}_no_stop_3.fastq

echo -e "[SGE - $(date +"%T")]\tBuilding .dict file"
java -jar /share/ClusterShare/software/contrib/briglo/picard/build/libs/picard.jar CreateSequenceDictionary R=${REF2}

grep stop_gained ${STEM}_variants_3.csq.vcf

for POS in $(grep stop_gained ${STEM}_variants_3.csq.vcf | awk '{print $2}'); do BASE=$(grep ${POS} ${STEM}_variants_3.csq.vcf | awk '{print $5}'); samtools view -h ${BAM3} | sam2tsv.jar -r ${REF2} | awk -v a=${POS} -v b=${BASE} '($8==a) && ($6==b)'; done > ${STEM}_reads_to_filter_3.tsv
python3 ${WORK_DIR}/scripts/filter_stops.py -f ${FASTQ3} -t ${STEM}_reads_to_filter_3.tsv -w ${STEM}_with_stop_3.fastq -n ${STEM}_no_stop_3.fastq
POS=''
BASE=''
# ${WORK_DIR}/scripts/filter_stops.sh ${WORK_DIR} 3 ${STEM} ${REF2} ${BAM3} ${STEM}_variants_3.csq.vcf ${FASTQ3}

# map fastq files and index bam files


FASTQ4=${STEM}_no_stop_3.fastq

# minimap2 -x map-ont -t 8 -k15 ${REF2} ${FASTQS}/${BASE%.srt*}.fastq > ${STEM}_no_stop_2.paf
# minimap2 -ax map-ont -k15 -t 8 ${REF2} ${FASTQS}/${BASE%.srt*}.fastq | samtools view -Sb - | samtools sort -o ${STEM}_no_stop_2.srt.bam -

minimap2 -cx map-ont -t 8 -k15 ${REF2} ${FASTQ4} > ${STEM}_pt_ref_no_stop_mapped.paf
minimap2 -ax map-ont -k15 -t 8 ${REF2} ${FASTQ4} | samtools view -Sb - | samtools sort -o ${STEM}_pt_ref_no_stop_mapped.srt.bam -

BAM4=${STEM}_pt_ref_no_stop_mapped.srt.bam
samtools index ${BAM4}


echo -e "[SGE - $(date +"%T")]\tDONE"







# R commands to get deletion site mapping points

#
# NanoFilt -q 10
