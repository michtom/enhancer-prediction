import numpy as np
import pandas as pd
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
from itertools import product


def read_bed(filename):
    content = []
    with open(filename) as f:
        for line in f:
            content.append(line.strip().split())
    return np.array(content)

def get_chr_names(start, stop):
    chromosomes_names = []
    prefix = "chr"
    for i in range(start, stop+1):
        chromosomes_names.append(prefix + f"{i}")
    return chromosomes_names

def get_chr_index_dict(DNA_array, chr_names):    
    chromosomes_dict = {x: 0 for x in chr_names}
    for i in range(len(DNA_array)):
        if DNA_array[i].id in chromosomes_dict:
            chromosomes_dict[DNA_array[i].id] = i
        else:
            continue
    return chromosomes_dict

def is_n_percent_greater_5(seq):
    if seq.count("N")/len(seq) > 0.05:
        return True
    else:
        return False

def split_chromosome(chr_name, chr_dict, DNA_array, enh_len):
    index = chr_dict[chr_name]
    seq = DNA_array[index].seq
    chr_df = pd.DataFrame()
    n = len(seq)
    for i in tqdm(range((n//enh_len)+1)):
        start = i*enh_len
        stop = min((i+1)*enh_len, n) - 1
        seq_part = str(seq[start:stop+1])
        if is_n_percent_greater_5(seq_part):
            data = pd.DataFrame({"start": start, "stop": stop, "seq": ""}, index=[0])
            chr_df = chr_df.append(data, ignore_index=True)
        else:
            seq_part = seq_part.replace("N", "")
            data = pd.DataFrame({"start": start, "stop": stop, "seq": seq_part}, index=[0])
            chr_df = chr_df.append(data, ignore_index=True)

    filename =  chr_name + ".csv"
    chr_df.to_csv(filename, index=False)
    del chr_df

def inverse_mer(mer):
    inv_mer = ''
    n = len(mer)
    for i in range(n-1, -1, -1):
        if(mer[i]=='A'):
            inv_mer += 'T'
            continue
        elif(mer[i]=='T'):
            inv_mer += 'A'
            continue
        elif(mer[i]=='C'):
            inv_mer += 'G'
            continue
        elif(mer[i]=='G'):
            inv_mer += 'C'
            continue
    return inv_mer

def count_mers(sequence, k):
    counts = {''.join(key): 0 for key in list(product('ACTG', repeat = k))}
    for i in range(len(sequence)-k+1):
        counts[sequence[i:(i+k)]] += 1
    for key in counts.copy():
        inv_mer = inverse_mer(key)
        val = counts[inv_mer]
        if(inv_mer != key):
            try:
                counts[key] += val
                del counts[inv_mer]
            except:
                continue
        else:
            continue
    n = len(sequence)
    for key in counts:
        counts[key] = counts[key]/n if n!= 0 else n
    
    return counts

def get_mers_dataframe(df, chr_name, k):
    kmers_df = pd.DataFrame()
    for i in tqdm(range(df.shape[0])):
        data = {"start": df.iloc[i,0], "stop": df.iloc[i,1]}
        seq = df.iloc[i,2]
        if type(seq) == type("abc"):
            k_mers = count_mers(seq, k)
        else:
            k_mers = count_mers("", k)
        data.update(k_mers)
        kmers_df = kmers_df.append(data, ignore_index=True)
    
    filename = f"{k}" + "_mers_" + chr_name + ".csv"
    kmers_df.to_csv(filename, index=False)
    del kmers_df