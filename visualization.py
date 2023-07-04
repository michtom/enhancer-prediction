import os

import plotly.express as px
import pandas as pd
from plotly.offline import plot
import xgboost as xb
import numpy as np


def calculate_predictions(chrom, model):
    # calculate percentage preditions for each window in chromosome (and replace unknown with mean)
    indices_unknown = chrom.index[chrom.sum(axis=1) == 0].tolist()
    df_without_unknown = chrom[~chrom.index.isin(indices_unknown)]
    avg_prob = 100 * model.predict_proba(df_without_unknown)[:, 1].mean()
    predictions = 100 * model.predict_proba(chrom)[:, 1]
    predictions[indices_unknown] = avg_prob
    return predictions


def make_visualization(probabilities, percentage, start_index=0):
    windows = range(1, len(probabilities))
    x_axis = windows[start_index:start_index+int((percentage/100)*len(probabilities))]
    y_axis = probabilities[start_index:start_index+int((percentage/100)*len(probabilities))]
    color = np.where(y_axis>50, 'possibly enhancer', 'possibly not an enhancer')
    df = pd.DataFrame({"x_axis": x_axis, "y_axis": y_axis, "color": color})
    fig = px.bar(df, x="x_axis", y="y_axis",
                 title="Chance of being an enhancer, in %",
                 width=1500, height=600,
                 color='color',
                 color_discrete_map={"possibly enhancer":"green", "possibly not an enhancer": "red"},
                 labels={"x_axis": "Number of window in chromosome", "y_axis": "percent of chance",
                         "color": "prediction"})
    fig.update_layout(title_x=0.5)
    fig.update_layout(showlegend=False)
    return plot(fig, output_type='div')


def get_visualization(k_mers_fpath, model, percent=1.0, start_index=0, predict_only=False):

    df = pd.read_csv(k_mers_fpath)
    X = df.drop(columns=["start", "stop"])
    predictions = calculate_predictions(X, model)
    df = pd.DataFrame([predictions])
    path = os.path.abspath(__file__).split("\\")
    path.pop()
    df.to_csv(os.path.join(os.getcwd(), "static", "tmp", "results.csv"), index=False, header=False)
    df.to_csv("results.csv", index=False, header=False)
    if not predict_only:
        plot = make_visualization(predictions, percent, start_index)
        return plot

