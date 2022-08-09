import pandas as pd
import numpy as np
import os
from joblib import Parallel, delayed
import shutil
import sys


def hier_dtype(path):
    exp_rows = np.array(pd.read_csv(path, index_col=0, nrows=0).columns.tolist())
    subj_ids = [x for x in exp_rows if "GTEX" in x]
    col_types = {"CHR": str, "POS": int, "VARIANT_ID": str, "GENE_ID": str, "TISSUE_ID":str, "VARIANT_ANNOTATION": str}

    for subj_id in subj_ids:
        col_types[subj_id] = str
    
    return col_types


def hier_join(i, ct_start, pair_chunks, hier_level):
    if hier_level == 0:
        path = base_path
    if hier_level > 0:
        path = hier_path
    counter = ct_start
    for pair in pair_chunks[i]:
        f1, f2 = pair
        
        ctypes1 = hier_dtype(f"{path}{f1}")
        ctypes2 = hier_dtype(f"{path}{f2}")

        df1 = pd.read_csv(f"{path}{f1}", dtype=ctypes1)
        df2 = pd.read_csv(f"{path}{f2}", dtype=ctypes2)

        # merge dataframes
        mdf = df1.merge(df2, on=grp_vars, how="outer", left_index=False, right_index=False)
        # write out merged dataframe with hierarchy and counter as denomination
        mdf.to_csv(f"{hier_path}{counter}.{hier_level}.csv", index=False)
        counter += 1


# base_path: path to pivot csv files from GeneSubjectCol.py
# hier_path: path to hold intermediate hierarchically joined data
def run_join(base_path, hier_path, hier_start=0):
    # make the directory if it does not exist yet
    if not os.path.exists(hier_path):
        os.makedirs(hier_path)

    hier_level = 0
    next_iter_len = -1
    while next_iter_len != 1:
        # hier_level of 0 will be the first files to merge
        # hier_level greater than 0 will have their own directory
        if hier_level == 0:
            file_paths = [x for x in os.listdir(base_path) if "GTEX" in x]
            # If there are an odd number of files to pair, defer the left-out file to next hierarchy
            if len(file_paths) % 2 != 0:
                shutil.copy(f"{base_path}{file_paths[-1]}", f"{hier_path}x0.{hier_level}.csv")

        if hier_level > 0:
            file_paths = [x for x in os.listdir(hier_path) if x.split(".")[1] == str(hier_level - 1)]
            if len(file_paths) % 2 != 0:
                shutil.copy(f"{hier_path}{file_paths[-1]}", f"{hier_path}x0.{hier_level}.csv")

        if len(file_paths) == 1:
            return hier_level

        # pair files with each other to merge
        pair_list = list(zip(*[iter(file_paths)]*2))
        # chunk pairs into no more than 10 chunks for parallel processing
        n_chunks = min((len(pair_list) // 2) + len(pair_list) % 2, 10)
        pair_chunks = np.array_split(pair_list, n_chunks)

        # parallel cannot access a global counter variable. ct_starts supplements this by providing counter start values
        ct_starts = np.array([len(x) for x in pair_chunks])
        ct_starts = np.insert(ct_starts, 0, 0)[:-1].cumsum()

        Parallel(n_jobs=-1)(delayed(hier_join)(i, ct_start, pair_chunks, hier_level) for i, ct_start in zip(range(n_chunks), ct_starts))

        hier_level += 1
        

if __name__ == "__main__":
    out_name = sys.argv[1]
    fname = os.path.split(out_name)[1].split(".")[0]
    # Declare important paths
    hier_path = f"./data/hier_merge_{fname}/"
    base_path = f"./data/pivots_subset_{fname}/"
    
    # universal columns to group by
    grp_vars = ["CHR", "POS", "VARIANT_ID", "GENE_ID", "TISSUE_ID", "VARIANT_ANNOTATION"]

    # Run hierarchical join and record final hierarchy level
    final_hier = run_join(base_path, hier_path)

    # Move final file to present directory
    shutil.move(f"{hier_path}x0.{final_hier}.csv", f"{out_name}")

    # cleanup after self
    shutil.rmtree(hier_path)
    shutil.rmtree(base_path)
