import pandas as pd


def subset_sum(ss_cols, ind):
    for_colns = pd.read_csv("./transcript_tpm.gct", sep = "\t", skiprows=2, usecols=ss_cols, low_memory=True)
    
    out_df = pd.concat([for_colns.iloc[:, :2], (for_colns[ss_cols[2:]] >= 1).sum(axis=1)], axis=1).rename(columns={0: "count"})
    
    out_df.to_csv(f"./intermed/{ind}.csv", index=False)

    return


subj_ids = pd.read_csv("./transcript_tpm.gct", sep = "\t", skiprows=2, nrows=0).columns

ss_cols_list = []
for n in range((len(subj_ids) // 2000) + 1):
    ss_cols_list.append(["transcript_id", "gene_id"] + list(subj_ids[2 + (n * 2000): min(((n * 2000) + 2 + 2000), len(subj_ids))]))

for i, ss_cols in enumerate(ss_cols_list):
    subset_sum(ss_cols, i)