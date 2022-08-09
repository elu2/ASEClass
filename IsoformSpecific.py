import pandas as pd
import sys
import os


def single_isoform_variants(var_df, tissue=None, genes=None):
    isoform_spec = var_df[var_df.transcript_count == 1]
    if tissue != None:
        isoform_spec = isoform_spec[isoform_spec.TISSUE_ID.isin(tissue)]
    
    if genes != None:
        isoform_spec = isoform_spec[isoform_spec.GENE_ID.isin(genes)]

    isoform_spec = isoform_spec.dropna(how="all", axis=1)
    
    return isoform_spec


if __name__ == "__main__":
    in_name = sys.argv[1]
    file_name = os.path.split(in_name)[1]
    
    var_df = pd.read_csv(in_name, sep="\t")

    single_isoform_variants(var_df).to_csv(f"{file_name.split('.')[0]}.isoform_specific.tsv",
                                           sep="\t", index=False)
