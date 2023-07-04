from django import forms
from django.core.validators import MaxValueValidator, MinValueValidator


class InputGeneObjectForm(forms.Form):
    chromosome_num = forms.IntegerField()

    def __init__(self, *args, **kwargs):
        super(InputGeneObjectForm, self).__init__(*args, **kwargs)
        self.fields['chromosome_num'].label = ""


class PlotPropertiesForm(forms.Form):
    chromosome_num = forms.IntegerField()
    min_index = forms.IntegerField()
    percentage = forms.FloatField()

    def __init__(self, *args, **kwargs):
        super(PlotPropertiesForm, self).__init__(*args, **kwargs)
        self.fields['min_index'].label = "Choose the first window index: "
        self.fields['chromosome_num'].label = "Choose the chromosome: "
        self.fields['percentage'].label = "Choose the percentage of data to visualize: "

