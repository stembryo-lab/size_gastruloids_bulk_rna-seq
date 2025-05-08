# Gene Markers UMAPs
In this folder are founf the scripts for downloading, processing and generating marker UMAPs from public embryo and gastruloid scRNA-Seq datasets.

- [Gastruloid dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229513)
- [Embryo dataset](https://marionilab.github.io/ExtendedMouseAtlas/)

## 0. Obtain data
To download and store datasets, run:

`./download_datasets.sh`

## 1. Setup environemnt

We advise to ue conda to create a separate environemtn with all the required packages to reproduce the results.
You can setup the environment as:

```bash
conda env create -n NAME_OF_ENVIRONMENT -f environment.yml
```

and activate it using 

```
conda activate NAME_OF_ENVIRONMENT
```

## 2. Analysis

Run `analysis_markers_umaps_ANALYSIS.ipynb` notebook to preprocess, analyze and generate UMAPs for markers in *markers.txt*.
