# Creates exon_region_classification

import pandas as pd
import numpy as np
import os
import itertools
from joblib import Parallel, delayed
from operator import itemgetter


def df_chunker(full_df, chunks):
    dfs = list()
    interval_size = full_df.shape[0]//chunks
    dfs.append(full_df.iloc[0:interval_size, :])

    for i in range(chunks - 1):
        dfs.append(full_df.iloc[(interval_size * (i + 1))
                   :(interval_size * (i + 2)), :])

    if full_df.shape[0] % chunks != 0:
        dfs.append(full_df.iloc[interval_size * chunks:, :])

    return dfs


def inter_ranges(intersect_counts, start_bp):
    counts = []
    cons_counts = [0]
    for count, it in itertools.groupby(intersect_counts):
        counts.append(count)
        cons_counts.append(len(list(it)))

    bp_is = np.cumsum(cons_counts) + start_bp

    ranges = np.lib.stride_tricks.sliding_window_view(bp_is, 2).copy()
    ranges[:, 1] = ranges[:, 1] - 1
    
    return counts, ranges


def get_transc_id(exon_df, bp_range):
    rstart, rend = bp_range
    ol_df = exon_df[(exon_df["start"] <= rstart) & (exon_df["end"] >= rend)]
    transc_ids = ";".join(ol_df["transcript_id"])
    return transc_ids
        

def range_classifier(model_exon_chunk):
    for i in range(model_exon_chunk.shape[0]):
        model_row = model_exon_chunk.iloc[i]
        # Get model data
        gene_id = model_row["gene_id"]
        chrom = model_row["seqname"]

        # Standardize ranges to be indices as well
        model_exon_range = np.arange(model_row["start"], model_row["end"] + 1) - model_row["start"]

        # Limit searching ranges to just same gene
        gene_exons = pri_exons[pri_exons["gene_id"] == gene_id]
        max_transcs = len(set(gene_exons.transcript_id))

        # Maintain number of isoforms overlapping with model exons' regions
        intersect_counts = np.zeros(len(model_exon_range))
        for j in range(gene_exons.shape[0]):
            pri_row = gene_exons.iloc[j]
            pri_exon_range = np.arange(pri_row["start"], pri_row["end"] + 1) - model_row["start"]
            intersection = np.intersect1d(model_exon_range, pri_exon_range)

            intersect_counts[intersection] = intersect_counts[intersection] + 1

        counts, ranges = inter_ranges(intersect_counts, model_row["start"])
        t_ids = [get_transc_id(gene_exons, x) for x in ranges]

        starts = list(map(itemgetter(0), ranges))
        ends = list(map(itemgetter(1), ranges))

        max_counts = [max_transcs] * len(counts)
        gene_ids = [model_row["gene_id"]] * len(counts)
        chroms = [model_row["seqname"]] * len(counts)

        dat_df = pd.DataFrame({
            "seqname":chroms,
            "start": starts,
            "end": ends,
            "transcript_count": counts,
            "max_transcripts": max_counts,
            "gene_id": gene_ids,
            "transcript_ids": t_ids}
        )

        dat_df.to_csv(
            "./data/exon_region_classification.csv", mode='a', index=False,
            header=not os.path.exists("./data/exon_region_classification.csv")
        )


if __name__ == "__main__":
    model_exons = pd.read_csv("./data/pcModelExons.csv")
    pri_exons = pd.read_csv("./data/pcPriExons.csv")
    
    chunks = df_chunker(model_exons, 100)

    par_out = Parallel(n_jobs=-1)(delayed(range_classifier)(mec) for mec in chunks)
