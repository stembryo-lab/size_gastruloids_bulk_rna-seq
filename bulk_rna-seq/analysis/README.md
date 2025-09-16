# Downstream Analysis

Run the `sizeproject_analysis.ipynb` notebook to perform the complete downstream analysis described in the publication and to reproduce the associated figures.

> **Note:** Make sure you have the Python modules listed in `requirements.txt` (or the Docker image) installed in your environment beforehand.

---

## Repository Files

| File | Description |
|------|-------------|
| **custom_functions.py** | Python module containing helper functions used across notebooks and scripts (data loading, plotting utilities, QC summaries, etc.). |
| **markers.xlsx** | Excel sheet listing marker genes or gene sets of interest for downstream analysis or validation. |
| **parameters.json** | JSON configuration file storing pipeline parameters (settings, thresholds, and paths). |
| **sizeproject_analysis.ipynb** | Jupyter notebook performing data exploration and downstream analysis on the processed RNA-seq count data. |
