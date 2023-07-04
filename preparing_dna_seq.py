import numpy as np
import pandas as pd
from Bio import SeqIO
import gzip
import warnings

warnings.filterwarnings('ignore')
from preprocessing_functions import *


def main():
    positive = read_bed("GM12878.bed")
    positive_df = pd.DataFrame(positive)
    positive_df = positive_df.loc[~positive_df[0].isin(["chrX"])]

    positive_df[1] = pd.to_numeric(positive_df[1])
    positive_df[2] = pd.to_numeric(positive_df[2])

    positive_df["length"] = positive_df[2] - positive_df[1] + 1
    enh_len = round(np.mean(positive_df.iloc[:, 4]), -2)

    DNA_array = []
    with gzip.open("GRCh37.primary_assembly.genome.fa.gz", "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            DNA_array.append(record)

    chromosomes_names = get_chr_names(1, 22)
    chromosomes_dict = get_chr_index_dict(DNA_array, chromosomes_names)

    for chr_name in chromosomes_names:
        split_chromosome(chr_name, chromosomes_dict, DNA_array, enh_len)


main()
