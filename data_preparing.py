from collections import Counter
import numpy as np
import pandas as pd
import Bio
from Bio import SeqIO
import gzip
import xgboost as xgb
from itertools import product
from tqdm import tqdm_notebook as tqdm
from sklearn.metrics import accuracy_score, precision_score, roc_auc_score, f1_score, confusion_matrix, recall_score


def inverse_mer(mer):
    # change a mer to its inverse version
    inv_mer = ''
    n = len(mer)
    for i in range(n - 1, -1, -1):
        if mer[i] == 'A':
            inv_mer += 'T'
            continue
        elif mer[i] == 'T':
            inv_mer += 'A'
            continue
        elif mer[i] == 'C':
            inv_mer += 'G'
            continue
        elif mer[i] == 'G':
            inv_mer += 'C'
            continue
    return inv_mer


def count_mers(sequence, k):
    # count how many times each mer of length k occured
    counts = {''.join(key): 0 for key in list(product('ACTG', repeat=k))}
    for i in range(len(sequence) - k + 1):
        counts[sequence[i:(i + k)]] += 1
    for key in counts.copy():
        inv_mer = inverse_mer(key)
        val = counts[inv_mer]
        if (inv_mer != key):
            try:
                counts[key] += val
                del counts[inv_mer]
            except:
                continue
        else:
            continue
    n = len(sequence)
    for key in counts:
        counts[key] = counts[key] / n if n != 0 else n

    return counts


def read_bed(filename):
    # read .bed files
    content = []
    with open(filename) as f:
        for line in f:
            content.append(line.strip().split())
    return np.array(content)


def get_mers_dataset(dataset, k, positive=True):
    # return a dataset with counted k-mers occurrence
    final_df = pd.DataFrame()
    for seq in tqdm(dataset["sequence"][:10]):
        seq = seq.replace("N", "")
        res_dict = count_mers(seq, k)
        if positive:
            res_dict["enhancer"] = 1
        else:
            res_dict["enhancer"] = 0
        final_df = final_df.append(res_dict, ignore_index=True)

    return final_df.fillna(0)


def read_fasta(file):
    # read .fasta files
    array = []
    with gzip.open(file, "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            array.append(record)


def get_train_test_datasets(positive_bed, negative_bed, array):
    """
    :param positive_bed: dataset with positive sequences from .bed file
    :param negative_bed: dataset with negative sequences from .bed file
    :param array: array of records from fasta file (SeqRecord objects)
    :return: tuple of 3 created datasets: with positive sequences for training, negative sequences for training
    and positive sequences for testing
    """
    positive_df = pd.DataFrame(positive_bed)
    negative_df = pd.DataFrame(negative_bed)
    positive_df = positive_df.loc[~positive_df[0].isin(["chrX"])]
    negative_df = negative_df.loc[~negative_df[0].isin(["chrX", "chrY"])]

    positive_df = positive_df.loc[~positive_df[0].isin(["chrX"])]
    positive_test = positive_df.loc[positive_df[0].isin(["chr1", "chr14", "chr21"])]
    positive_train = positive_df.loc[~positive_df[0].isin(["chr1", "chr14", "chr21"])]
    negative_train = negative_df.loc[~negative_df[0].isin(["chrX", "chrY"])]

    array_dict = {}
    for i in range(len(array)):
        array_dict[array[i].name] = array[i].seq
    newcol = []
    for row in positive_train.iterrows():
        name = row[1][0]
        start = row[1][1]
        stop = row[1][2]
        frag = array_dict[name][int(start) - 1:int(stop)]
        newcol.append(frag)

    positive_train["sequence"] = newcol

    newcol = []
    for row in positive_test.iterrows():
        name = row[1][0]
        start = row[1][1]
        stop = row[1][2]
        frag = array_dict[name][int(start) - 1:int(stop)]
        newcol.append(frag)

    positive_test["sequence"] = newcol

    newcol = []
    for row in negative_train.iterrows():
        name = row[1][0]
        start = row[1][1]
        stop = row[1][2]
        frag = array_dict[name][int(start) - 1:int(stop)]
        newcol.append(frag)

    negative_train["sequence"] = newcol
    return positive_train, negative_train, positive_test


def get_avg_enhancer_length(positive_train, positive_test):
    # return an average length of enhancer, in both train and test datasets, rounded to 100
    all_positives = pd.concat([positive_train, positive_test])
    avg_enh_length = all_positives["sequence"].str.len().mean()
    return round(avg_enh_length, -2)
