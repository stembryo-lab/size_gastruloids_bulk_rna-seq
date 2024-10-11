#!/bin/bash
# Download raw fastq files and transcriptome references

# create folders
FASTQS_DIR='data/fastqs/raw'
REF_DIR='data/references'

# urls
FASTA_URL='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz'
ANNOT_URL='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz'

if [[ ! -d $FASTQS_DIR ]]
then
    mkdir -p data/fastqs/raw
fi

if [[ ! -d $REF_DIR ]]
then
    mkdir -p data/references/{fasta,annotation}
fi

# download files
# fastqs
#wget [link2fastqs] -o $FASTQS_DIR

# references
wget $FASTA_URL -P $REF_DIR/fasta
wget $ANNOT_URL -P $REF_DIR/annotation