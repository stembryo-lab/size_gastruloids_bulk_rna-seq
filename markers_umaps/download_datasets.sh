#!/bin/bash
# Script that downloads scRNA-Seq datasets

# create download folders
mkdir -p data/{mouse,gastruloid}

# URLs
# Mouse: https://journals.biologists.com/dev/article/151/3/dev201867/342647/Tracking-early-mammalian-organogenesis-prediction
PIJUANSALA_MOUSE_URL='https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz'
# Gastruloid: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229513
PRISCA_GASTRULOID_URL=https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE229513&format=file&file=GSE229513%5Fgastruloidsobject%2Erds%2Egz
# PRISCA_GASTRULOID_URL1='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE229513&format=file&file=GSE229513%5FUMI%5Fcounts%2Emtx%2Egz'
# PRISCA_GASTRULOID_URL2='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE229513&format=file&file=GSE229513%5Fbarcodes%2Etsv%2Egz'
# PRISCA_GASTRULOID_URL3='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE229513&format=file&file=GSE229513%5Fgenes%2Etsv%2Egz'

mkdir -p data/mouse
wget $PIJUANSALA_MOUSE_URL -P data/mouse/
cat data/mouse/atlas_data.tar.gz data/mouse/atlas_data.tar.gz.1 > data/mouse/atlas_data.tar.gz
gunzip data/mouse/atlas_data.tar.gz
tar -xvf data/mouse/atlas_data.tar -C data/mouse
cp -r data/mouse/atlas/ data/mouse/
rm -rf data/mouse/atlas
rm -rf data/mouse/*.tar*

mkdir -p data/gastruloid
wget $PIJUANSALA_MOUSE_URL -O data/gastruloid/GSE229513_gastruloidsobject.rds.gz
# wget $PRISCA_GASTRULOID_URL1 -O data/gastruloid/GSE229513_UMI_counts.mtx.gz
# wget $PRISCA_GASTRULOID_URL2 -O data/gastruloid/GSE229513_barcodes.tsv.gz
# wget $PRISCA_GASTRULOID_URL3 -O data/gastruloid/GSE229513_genes.tsv.gz
gunzip data/gastruloid/GSE229513_gastruloidsobject.rds.gz
# gunzip data/gastruloid/GSE229513_UMI_counts.mtx.gz
# gunzip data/gastruloid/GSE229513_barcodes.tsv.gz
# gunzip data/gastruloid/GSE229513_genes.tsv.gz
