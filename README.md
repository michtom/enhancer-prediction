# wb_2023

A django web service that predicts and visualises the probability that a DNA sequence is an enhancer. The user can choose a chromosome and a range of solutions displayed. Predictions made with xgboost, visualisations created using plotly.

Project made during Case Study 2 course at Faculty of Mathematics and Information Science, Warsaw University of Technology.

Authors:
* Jędrzej Sokołowski
* Filip Szympliński
* Michał Tomczyk

Functions in `data_preparing.py`, `preparing_dna_kmers.py`,  `preparing_dna_seq.py` and `preprocessing_functions.py` were run on raw .bed and .fa files. Preprocessed data, divided into .csv files for each chromosome can be downloaded from https://drive.google.com/drive/folders/1xTG9NMLiQNMETCdwUvtM6oJZ71c7mMuh. To make the website run properly, fill enhancer-prediction/wb_site/static/4_mers/ folder with this unpacked files.

