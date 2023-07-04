import pandas as pd
import warnings
warnings.filterwarnings('ignore')
from threading import Thread
import psutil
import time
from preprocessing_functions import *

def threds_working(chromosomes_names):
    threads = []
    for chr_name in chromosomes_names:
        filename = chr_name + ".csv"
        df = pd.read_csv(filename, sep=",")
        x = Thread(target=get_mers_dataframe, args=(df, chr_name, 4,), )
        threads.append(x)
        x.start()
        while psutil.virtual_memory().percent > 79:
            print("WĄTEK CZEKA: ", x)
            time.sleep(5)

    for t in threads:
        t.join()
    
    print("All done!")

def main():
    chromosomes_names = get_chr_names(1, 4)
    threds_working(chromosomes_names)

main()