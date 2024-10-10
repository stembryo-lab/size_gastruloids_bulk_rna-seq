#!/bin/bash

mkdir -p data/{mouse,gastruloid}

# Download scRNA-Seq datasets

# Mouse: https://journals.biologists.com/dev/article/151/3/dev201867/342647/Tracking-early-mammalian-organogenesis-prediction
PIJUANSALA_MOUSE_URL='https://cloud.mrc-lmb.cam.ac.uk/s/yxq7FRtYsLyF3jQ/download'
# Gastruloid: https://w'ww.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229513
PRISCA_GASTRULOID_URL='https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE229513&format=file'

wget $PIJUANSALA_MOUSE_URL -P data/mouse/
unzip data/mouse/download
rm -rf download
mv data/mouse/ExtendedMouseAtlas/embryo_complete.h5ad ../data/mouse
rm -rf data/mouse/ExtendedMouseAtlas # keep just .h5ad file

wget $PRISCA_GASTRULOID_URL

