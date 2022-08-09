import pandas as pd
import os

grp_vars = ["CHR", "POS", "VARIANT_ID", "GENE_ID", "TISSUE_ID", "VARIANT_ANNOTATION"]

f_path = "./data/pivots/"
fnames = [x for x in os.listdir(f_path) if "GTEX" in x]

base_df = pd.read_csv(f"{f_path}{fnames[0]}")
counter = 0
for f in fnames[1:]:
    to_merge = pd.read_csv(f"{f_path}{f}")
    mdf = base_df.merge(to_merge, on=grp_vars, how="outer", left_index=False, right_index=False)
    base_df = mdf
    print(counter)
    counter += 1

base_df.to_csv("./data/full_merged.csv", index=False)