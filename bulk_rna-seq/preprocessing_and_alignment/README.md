# Bulk RNA-seq Preprocessing, QC and Mapping

This repository contains a simple, reproducible pipeline to process bulk RNA-seq data from raw FASTQ files through quality control, adapter/quality trimming, and alignment to a reference genome.

---

## Pipeline Overview

1. **Data download**  
   `download_data.sh` retrieves raw FASTQ files from public repositories (e.g. SRA/ENA) or institutional storage.

2. **Quality control & trimming**  
   `qc_and_trim.sh` runs initial quality checks with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and performs adapter/quality trimming with [fastp](https://github.com/OpenGene/fastp) or [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).  
   - Generates pre- and post-trim QC reports.  
   - Ensures only high-quality reads proceed to mapping.

3. **Alignment (mapping)**  
   `mapping.sh` aligns trimmed reads to a reference genome using a splice-aware aligner such as [STAR](https://github.com/alexdobin/STAR) or [HISAT2](https://daehwankimlab.github.io/hisat2/).  
   - Produces raw countdata per sample (ReadsPerGene.out.tab)
  
4. **Concatenate countdata**
   
   `create_countdata.py` concatenates counts from each aligned sample and creates a merged count matrix to be used in the downstream analysis part.

