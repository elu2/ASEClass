import pandas as pd
import numpy as np
import os
from joblib import Parallel, delayed
import sys


def in_region(qpos, qseqname):
    val_betw = np.where(
        ((qpos - class_df["start"]) * (class_df["end"] - qpos)) >= 0)
    val_df = class_df.iloc[val_betw]
    try:
        found = val_df[val_df["seqname"] == qseqname].iloc[0]
    except IndexError:
        return None, None, None, None

    return found["transcript_count"], found["max_transcripts_exon"], found["max_transcripts_gene"], found["transcript_ids"]


def variant_classifier(sdf):
    for i in range(sdf.shape[0]):
        samp_row = sdf.iloc[i]

        qseqname = samp_row["CHR"]
        qpos = samp_row["POS"]
        qVID = samp_row["VARIANT_ID"]

        tc, mte, mtg, tid = in_region(qpos, qseqname)

        row = pd.DataFrame({"VARIANT_ID": [qVID], "transcript_count": [tc], "max_transcripts_exon": [mte],
                            "max_transcripts_gene": [mtg], "transcript_ids": [tid]})
        row.to_csv(f"{data_path}{fname_id}.variant_class_table.tsv", sep="\t", mode="a",
                   index=False, header=not os.path.exists(f"{data_path}{fname_id}.variant_class_table.tsv"))


def df_chunker(full_df, chunks):
    dfs = list()
    interval_size = full_df.shape[0]//chunks
    dfs.append(full_df.iloc[0:interval_size, :])

    for i in range(chunks - 1):
        dfs.append(full_df.iloc[(interval_size * (i + 1)): (interval_size * (i + 2)), :])

    if full_df.shape[0] % chunks != 0:
        dfs.append(full_df.iloc[interval_size * chunks:, :])

    return dfs


def post_run_join(ref_df, pr_df):
    mdf = ref_df.merge(pr_df, how="inner", on=["VARIANT_ID"])
    mdf = mdf.sort_values(["CHR", "POS"])
    return mdf


if __name__ == "__main__":
    # Read, process input path
    in_path = sys.argv[1]
    data_path, fname = os.path.split(in_path)
    data_path = data_path + "/"
    fname_id = fname.split(".")[0]

    # Read in data
    full_sample_df = pd.read_csv(f"{data_path}{fname}", sep=",")
    class_df = pd.read_csv("./data/exon_region_classification.csv")

    # Chunk data for parallel processing
    chunks = df_chunker(full_sample_df, 40)
    Parallel(n_jobs=-1)(delayed(variant_classifier)(sdf) for sdf in chunks)

    # Join reference and nacent
    pr_df = pd.read_csv(
        f"{data_path}{fname_id}.variant_class_table.tsv", sep="\t").drop_duplicates()

    # Merged dataframe
    mdf = post_run_join(full_sample_df, pr_df)
    # Overwrite nacent dataframe
    mdf.to_csv(f"{data_path}{fname_id}.variant_class_table.tsv",
               sep="\t", index=False)
