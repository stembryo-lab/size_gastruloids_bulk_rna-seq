# Gene Markers UMAPs
In this folder are founf the scripts for downloading, processing and generating marker UMAPs from public embryo and gastruloid scRNA-Seq datasets.

- [Gastruloid dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229513)
- [Embryo dataset](https://marionilab.github.io/ExtendedMouseAtlas/)

## 0. Obtain data
To download and store datasets, run:

`./download_datasets.sh`

## 1. Analysis
Run `analysis_markers_umaps.ipynb` notebook to preprocess, analyze and generate UMAPs for markers in *markers.txt*.

*Set:`dataset='gastruloid'` or `dataset='mouse'` to select which dataset to perform the analysis.*