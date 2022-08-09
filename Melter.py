import os
import pandas as pd


files = [x for x in os.listdir("./data/tissueSpecific/") if "variant_class_table" in x]

for file in files:
    tissue_id = file.split(".")[0]
    df = pd.read_csv(f"./data/tissueSpecific/{file}", sep="\t", dtype=str)
    df = df.dropna(subset="transcript_ids")

    df = df.melt(id_vars=["CHR", "POS", "GENE_ID", "VARIANT_ID", "TISSUE_ID", "VARIANT_ANNOTATION",
                     "transcript_count", "max_transcripts_exon", "max_transcripts_gene", "transcript_ids"],
            var_name="SUBJECT_ID",
            value_name="RATIO"
           )
    df = df.dropna(subset="RATIO")

    df.to_csv(f"./data/tissueSpecific/{tissue_id}.variant.melted.tsv", sep="\t", index=False)