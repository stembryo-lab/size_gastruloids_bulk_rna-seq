#!/bin/bash

#SBATCH -J star # A job-name
#SBATCH -p haswell # Partition 
#SBATCH -n 8 # Number of cores 
#SBATCH --mem 100000 #(GB per node)
#SBATCH -o ./jobs.out/slurm.%j.out
#SBATCH -e ./jobs.out/slurm.%j.err

: 'Mapping transcriptomic reads to genome with STAR'

#load modules in hpc
interactive
module load gzip/1.10-GCCcore-11.2.0 
module load STAR/2.7.10a-GCC-11.2.0 

# global options
DATASET_NAME=$1
DO_MAPPING=$2

# paths
PATH2READS='data/fastqs/raw'
PATH2FASTA='data/references/fasta/GRCm38.p6.genome.fa'
PATH2ANNOTATION='data/references/annotation/gencode.vM25.annotation.gtf'
PATH2GENOME_IDX='data/references/fasta/index'

MAPPING_OUT='star_out'
MULTIQC_OUT='star_out/reports'

# STAR parameters
THREADS=12
# index generation
GENOME_SA_INDEXBASES=12
GENOME_SA_SPARSE=3
#mapping
LOAD_MODE='NoSharedMemory'
READ_FORMAT='zcat'
GENE_ID='gene_name'
FILTER_MULTIMAP=10
OUTSAM_FORMAT='BAM Unsorted'
QUANTMODE='GeneCounts' 

cat << \
.
            MAPPING
=========================================

# PIPELINE GLOBAL OPTIONS
- Dataset Name: $DATASET_NAME
- Perform alignment: $DO_MAPPING

# DATA FOLDER
Raw fastq folder: $PATH2READS

# REFERENCES
Fasta: $(echo $PATH2FASTA | sed 's:.*/::')
Anntotation: $(echo $PATH2ANNOTATION | sed 's:.*/::')

# OUTPUT FOLDERS
Mapping output: $MAPPING_OUT
MultiQC report: $MULTIQC_OUT

# MAPPING PARAMETERS
# Memory
Threads: $THREADS
# memory for index generation
SA pre-indexingn length: $GENOME_SA_INDEXBASES
Indexing distance : $GENOME_SA_SPARSE
#mapping
Genome loading mode: $LOAD_MODE
Read format= $READ_FORMAT
Gene id in annotation: $GENE_ID
Multimappers filter: $FILTER_MULTIMAP
SAM/BAM Output format: $OUTSAM_FORMAT 
STAR quantMode: $QUANTMODE  


=========================================
.
\

if [ $# -ne 2 ]
then
    cat << \
.
    echo "Please, give:"
    echo "1) Naming for the dataset to analyze"   
    echo "2) "true/false" for performing Star Alignment"
.
\

fi


# 1. GENOME INDEX GENERATION
if [[ ! -f ${PATH2GENOME_IDX}/genomeParameters.txt ]] 
then
    mkdir -p ${PATH2GENOME_IDX}
    echo '
    Starting to index the genome...
    '

    STAR \
    --runMode genomeGenerate \
    --runThreadN $THREADS \
    --genomeDir $PATH2GENOME_IDX \
    --genomeFastaFiles $PATH2FASTA \
    --sjdbGTFfile $PATH2ANNOTATION \
    --genomeSAindexNbases $GENOME_SA_INDEXBASES \
    --genomeSAsparseD $GENOME_SA_SPARSE
else
    echo 'Genome index already computed.'
fi 


# 2.ALIGNMENT
if [[ $DO_MAPPING = true ]]
then
    for file in $PATH2READS/*_R1*.gz
    do
        sample=$(echo $file | sed 's:.*/::' | cut -d '_' -f 1)
        # single-end reads
        if [[ ! -f $PATH2READS/${sample}_R2*.gz ]]
        then
            r1=$PATH2READS/${sample}_R1*.fastq.gz
            if [[ ! -f $MAPPING_OUT/$sample/${sample}_Log.final.out ]]; then
            echo '
            Starting to map '${sample} 'reads
            '
                STAR \
                --genomeDir $PATH2GENOME_IDX \
                --genomeLoad $LOAD_MODE \
                --runThreadN $THREADS \
                --readFilesCommand $READ_FORMAT \
                --readFilesIn $r1  \
                --sjdbGTFfile $PATH2ANNOTATION \
                --sjdbGTFtagExonParentGene $GENE_ID \
                --outFilterMultimapNmax $FILTER_MULTIMAP \
                --outFileNamePrefix $MAPPING_OUT/$sample/${sample}_ \
                --outSAMtype $OUTSAM_FORMAT \
                --quantMode $QUANTMODE

            fi
            echo 'Mapping done for '${sample}

            # paired-end reads    
        elif [[ -f $PATH2READS/${sample}_R2*.gz ]]
        then
            r1=$PATH2READS/${sample}_R1*.gz
            r2=$PATH2READS/${sample}_R2*.gz
            if [[ ! -f $MAPPING_OUT/$sample/${sample}_Log.final.out ]]; then
                echo '
                Starting to map '${sample} 'reads
                '
                STAR \
                --genomeDir $PATH2GENOME_IDX \
                --genomeLoad $LOAD_MODE \
                --runThreadN $THREADS \
                --readFilesCommand $READ_FORMAT \
                --readFilesIn $r1 $r2 \
                --sjdbGTFfile $PATH2ANNOTATION \
                --sjdbGTFtagExonParentGene $GENE_ID \
                --outFilterMultimapNmax $FILTER_MULTIMAP \
                --outFileNamePrefix $MAPPING_OUT/$sample/${sample}_ \
                --outSAMtype $OUTSAM_FORMAT \
                --quantMode $QUANTMODE

            fi
            echo 'Mapping done for '${sample}
        fi

        # store mapping report in MultiQC folder
        mkdir -p ${MULTIQC_OUT}
        cp ${MAPPING_OUT}/${sample}/${sample}_ReadsPerGene.out.tab $MULTIQC_OUT/${sample}_ReadsPerGene.out.tab
    done
    module unload STAR/2.7.10a-GCC-11.2.0
    # run MultiQC on alignment reports
    module load MultiQC/1.9-foss-2019b-Python-3.7.4 
    multiqc $MULTIQC_OUT/* -n $DATASET_NAME -f -s -o $MULTIQC_OUT/

    echo '
    ==> STAR ALIGNMENT FINISHED <==
    --------------------------------

    Number of files analyzed: '$(ls ${PATH2READS}/*R1* -1 | wc -l)'
    Alignment output stored in: '${MAPPING_OUT}'
    MultiQC and STAR reports stored in: '${MULTIQC_OUT}'
    '
elif [[ $DO_MAPPING = false ]]
then
    echo 'No read mapping performed'
fi


