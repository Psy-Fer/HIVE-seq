
# HIVepsilon-seq - scalable characterisation of intact persistent proviral HIV reservoirs in women

## HIVe-seq

Run guppy basecaller to produce fastq file

If data is barcoded, dumultiplex and run on each barcode fastq

## Requirements:

python3.8 and java are required

## check java

```
java --version

# you should see something like the following:

openjdk version "11.0.16" 2022-07-19
OpenJDK Runtime Environment (build 11.0.16+8-post-Ubuntu-0ubuntu122.04)
OpenJDK 64-Bit Server VM (build 11.0.16+8-post-Ubuntu-0ubuntu122.04, mixed mode, sharing)
```

## Setup python environment

```
python3.8 -m venv HIVE-seq-py3.8
source HIVE-seq-py3.8/bin/activate


pip install --upgrade pip
pip install wheel setuptools cython
pip --no-cache-dir install pysamstats
```

## Running Hive-seq on basecalled fastq file

Argument order:

```
hive.sh ${FASTQ} ${WORK_DIR} ${REF} ${GFF} ${PICARD} > HIVE.log 2>&1
```

Run example on a single barcode fastq:

```
mkdir HIVE-analysis
cd HIVE-analysis

~/repos/HIVE-seq/hive.sh barcode11.fastq . ~/repos/HIVE-seq/ref/HXB2_trimmed.fasta ~/repos/HIVE-seq/ref/HXB2_AA_17Jul20.gff ~/install/picard-2.27.5/picard.jar > HIVE.log 2>&1

```
