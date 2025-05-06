#!/bin/bash
# Input variables
ANALYSIS_NAME=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
FASTQC_REPORT_DIR=$4

# parameters
TRIM_Q=20 # read quality threshold
TRIM_TYPE='--illumina'

cat << \
.
            QC & TRIMMING
=========================================

# PIPELINE GLOBAL OPTIONS
Analysis name: $ANALYSIS_NAME
# DATA FOLDER
Raw fastq folder: $INPUT_DIR

# OUTPUT FOLDERS
Trimmed reads output folder: $OUTPUT_DIR
FastQC and MultiQC reports folder: $FASTQC_REPORT_DIR

# TRIM GALORE PARAMETERS
Quality threshold: $TRIM_Q
Adapter type: $TRIM_TYPE

=========================================
.
\

if [ $# -ne 4 ]
then
    cat << \
.
    echo "Please, give:"
    echo "1) Naming for the dataset to analyze"   
    echo "2) Input data folder"   
    echo "3) Output data folder"
    echo "4) FastQC reports folder"
.
\

fi

: '
## QC ##

Run FastQC for every file in every condition. 
Then, MultiQC merges all FastQC reports.
'
# FastQC
mkdir -p ${FASTQC_REPORT_DIR}

echo '
PERFORMING FASTQC ANALYSIS
Starting FastQC analysis...
'
for file in ${INPUT_DIR}/R1.fastq.gz
do
    sample=$(echo $file | sed 's:.*/::' | cut -d '-' -f 1)

    for file in ${INPUT_DIR}/*.fastq.gz
    do
        fastqc -o ${FASTQC_REPORT_DIR} $file
        echo 'Analysis on '${sample}' done'
    done
done
echo '
==> FASTQC ANALYSIS FINISHED <==
--------------------------------

Number of files analyzed: '$(ls ${INPUT_DIR} -1 | wc -l)'
Reports stored in: '${FASTQC_REPORT_DIR}'

'
# MultiQC
echo 'MULTIQC
Merging FastQC reports...
'
multiqc ${FASTQC_REPORT_DIR}/* -n ${ANALYSIS_NAME}_raw -f -o $4 
echo '
==> MULTIQC ANALYSIS FINISHED <==
--------------------------------
Report name: '${ANALYSIS_NAME}'.html
Report stored in: '${FASTQC_REPORT_DIR} 

: '
## TRIMMING ##

Trim Illumina Adapter from read using TrimGalore!
'
echo '
PERFORMING ADAPTER TRIMMING
Running Trim Galore!...
'
for file in ${INPUT_DIR}/*R1.fastq.gz
do
    sample=$(echo $file | sed 's:.*/::' | cut -d '_' -f 1)

    # single-end samples 
    if [[ ! -f $INPUT_DIR/${sample}_R2.fastq.gz ]]
    then
        trim_galore --quality $TRIM_Q $TRIM_TYPE $INPUT_DIR/${sample}_R1.fastq.gz -o $OUTPUT_DIR
        echo 'Trimming for '${sample} ' done'

    elif [[ -f $INPUT_DIR/${sample}_R2.fastq.gz ]]
    then
        trim_galore --quality $TRIM_Q $TRIM_TYPE --paired $INPUT_DIR/${sample}_R1.fastq.gz $INPUT_DIR/${sample}_R2.fastq.gz -o $OUTPUT_DIR
        echo 'Trimming for '${sample} ' done'
    fi
done


    echo '
==> ADAPTER TRIMMING FINISHED <==
--------------------------------
Number of files analyzed: '$(ls ${INPUT_DIR} -1 | wc -l)'
Paired-end samples: '$(ls ${INPUT_DIR}/*R2* -1 | wc -l)'
Reports stored in: '${FASTQC_REPORT_DIR}'
'

echo '
PERFORMING FASTQC ANALYSIS
Starting FastQC analysis...
'
for file in ${OUTPUT_DIR}/*.fq.gz
do
    sample=$(echo $file | sed 's:.*/::' | cut -d '.' -f 1)

    fastqc -o $FASTQC_REPORT_DIR $file
    echo 'Analysis on '${sample}' done'

done
echo '
==> FASTQC ANALYSIS FINISHED <==
--------------------------------

Number of files analyzed: '$(ls ${INPUT_DIR} -1 | wc -l)'
Reports stored in: '${FASTQC_REPORT_DIR}'

'
# MultiQC
echo 'MULTIQC
Merging FastQC reports...
'
multiqc ${FASTQC_REPORT_DIR}/* -n ${ANALYSIS_NAME}_raw -f -o $FASTQC_REPORT_DIR 
echo '
==> MULTIQC ANALYSIS FINISHED <==
--------------------------------
Report name: '${ANALYSIS_NAME}'.html
Report stored in: '${FASTQC_REPORT_DIR} 


