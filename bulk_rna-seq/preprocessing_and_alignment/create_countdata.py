#!/usr/bin/env python3
import pandas as pd
import glob
import os

# Pattern to match your STAR count files
files = glob.glob("star_out/*/*ReadsPerGene.out.tab")  # adjust the path/pattern

# Pick which column of STAR to use:
# STAR columns: 0=GeneID, 1=unstranded, 2=1st-strand, 3=2nd-strand
column_index = 2  # 1-based counts column from STAR (2 = first-strand counts)

combined = None

for f in files:
    # Take sample name from filename (everything before first dot, adjust as needed)
    sample = os.path.basename(f).split('.')[0]
    # Read STAR file
    df = pd.read_csv(f, sep='\t', header=None)
    gene_ids = df.iloc[:, 0]
    counts = df.iloc[:, column_index]

    sample_df = pd.DataFrame({sample: counts})
    sample_df.insert(0, 'GeneID', gene_ids)

    if combined is None:
        combined = sample_df
    else:
        combined = combined.merge(sample_df, on='GeneID')

# Write combined matrix
combined.to_csv('sizeproject_counts.csv', sep=',', index=False)
print("Combined matrix written to sizeproject_counts.csv")
