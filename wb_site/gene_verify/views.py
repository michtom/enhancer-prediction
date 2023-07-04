import traceback

from django.http import HttpResponse
from django.shortcuts import render
from .forms import InputGeneObjectForm, PlotPropertiesForm
from django.utils.encoding import smart_str
from django.conf import settings
import os
import xgboost as xb
print(os.getcwd())
from visualization import get_visualization

# Create your views here.


def index(request):
    if request.method == "GET":
        return render(request, "index.html")


def home(request):
    if request.method == "GET":
        return render(request, "home.html")


def about(request):
    if request.method == "GET":
        return render(request, "about.html")


def classifier(request, action=''):

    if request.method == "POST":
        if 'button1' in request.POST:

            form = InputGeneObjectForm(request.POST)

            if form.is_valid():
                chromosome_num = form.cleaned_data["chromosome_num"]
                k_mers_fpath = os.path.join(settings.MERS_FOLDER_PATH, f"4_mers_chr{str(chromosome_num)}.csv")
                try:
                    xgb_model = xb.XGBClassifier()
                    xgb_model.load_model(os.path.join(settings.PROJ_DIR, "model.json"))
                    get_visualization(k_mers_fpath, xgb_model, predict_only=True) #, start_index=start_index)
                except Exception as e:
                    traceback.print_exc()

                with open(os.path.join(settings.TMP_FOLDER_PATH, "results.csv"), "r") as file:
                    data = file.read()

                response = HttpResponse(data, content_type='application/force-download')
                response['Content-Disposition'] = f'attachment; filename="results_chromosome{str(chromosome_num)}.csv"'

                return response

            return render(request, "classifier.html", {"form" : form, "predicted" : True})
        elif 'button2' in request.POST:
            form2 = PlotPropertiesForm(request.POST)
            if form2.is_valid():

                chromosome_num = form2.cleaned_data["chromosome_num"]
                percentage = form2.cleaned_data["percentage"]
                start_index = form2.cleaned_data["min_index"]
                k_mers_fpath = os.path.join(settings.MERS_FOLDER_PATH, f"4_mers_chr{str(chromosome_num)}.csv")
                try:
                    xgb_model = xb.XGBClassifier()
                    xgb_model.load_model(os.path.join(settings.PROJ_DIR, "model.json"))
                    data = get_visualization(k_mers_fpath, xgb_model, percentage, start_index=start_index)
                except Exception as e:
                    print(e)
                    data = 'error'
                return render(request, "classifier.html", {"form2": form2, 'plot_div': data, "plot": True})
        else:
            form = InputGeneObjectForm()
            form2 = PlotPropertiesForm()
            return render(request, "classifier.html", {"form": form, "form2" : form2})

    elif request.method == 'GET':
        if action == 'plot':
            form = InputGeneObjectForm(request.POST)

            form2 = PlotPropertiesForm()
            return render(request, "classifier.html", {"form": form, "form2": form2, "predicted": True, "plot": True})

    form = InputGeneObjectForm()
    form2 = PlotPropertiesForm()
    return render(request, "classifier.html", {"form" : form, "form2" : form2})



