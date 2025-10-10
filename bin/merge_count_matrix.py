#!/usr/bin/env python
# filepath: /home/mitch/Dokumente/Worklflow/Pipeline_Project/workflows_SS25_project/nf-core-workflowhunfeldruhland/bin/merge_count_matrix.py
import sys
import subprocess
import os
try:
    import pandas as pd
except ImportError:
    print("pandas not found, installing...", file=sys.stderr)
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas"])
    import pandas as pd



outdir = sys.argv[1]
tsv_files = sys.argv[2:]

merged_df = None

for f in tsv_files:
    df = pd.read_csv(f, sep="\t", comment='#')
    sample_col = [col for col in df.columns if col.endswith(".bam")][0]
    sample_name = sample_col.split(".")[0]
    
    if merged_df is None:
        merged_df = df.iloc[:, :1].copy()
        merged_df[sample_name] = df[sample_col]
    else:
        merged_df = pd.merge(merged_df, df[["Geneid", sample_col]], on="Geneid", how="outer")
        merged_df.rename(columns={sample_col: sample_name}, inplace=True)

os.makedirs(f"{outdir}/subread", exist_ok=True)
merged_df.to_csv(f"{outdir}/subread/merged_count_matrix.tsv", sep="\t", index=False)
print(f"✓ Merged count matrix saved to {outdir}/subread/merged_count_matrix.tsv")