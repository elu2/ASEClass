import pandas as pd
import numpy as np
import os
from joblib import Parallel, delayed
import sys
import shutil


# Configured for skeletal muscle for now
def pivot_samples(samp_files):
    for fname in samp_files:
        df = pd.read_csv(f"./GTEx_Analysis_v8_ASE_counts_by_subject/{fname}", sep="\t")
        df = df[df["TISSUE_ID"] == tissue]
        samp_id = fname.split(".")[0]
        samp_id_i = np.where(np.array(grp_vars + samp_ids) == samp_id)[0][0]

        df["CT_RATIO"] = df[["REF_COUNT", "ALT_COUNT"]].astype(str).agg("|".join, axis=1)
        sample_piv = df.pivot(index=grp_vars, columns="SUBJECT_ID", values="CT_RATIO")

        sample_piv.to_csv(f"./data/pivots_subset_{tissue}/{samp_id}.pivot.csv")


if __name__ == "__main__":
    tissue = sys.argv[1]
    
    if not os.path.exists(f"./data/pivots_subset_{tissue}/"):
        os.makedirs(f"./data/pivots_subset_{tissue}/")
    
    grp_vars = ["CHR", "POS", "VARIANT_ID", "GENE_ID", "TISSUE_ID", "VARIANT_ANNOTATION"]

    samp_files = [fname for fname in os.listdir("./GTEx_Analysis_v8_ASE_counts_by_subject") if "ase_table" in fname]
    samp_ids = [x.split(".")[0] for x in samp_files]
    samp_chunks = np.array_split(samp_files, 40)
    
    # Parallelize on samples
    Parallel(n_jobs=-1)(delayed(pivot_samples)(files_chunk) for files_chunk in samp_chunks)