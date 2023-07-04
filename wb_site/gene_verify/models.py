from django.db import models
from django.core.validators import MaxValueValidator, MinValueValidator
from decimal import Decimal

PERCENTAGE_VALIDATOR = [MinValueValidator(0), MaxValueValidator(100)]


# Create your models here.
class InputGeneObject(models.Model):
    chromosome_num = models.IntegerField(validators=[
        MaxValueValidator(22),
        MinValueValidator(1)
    ])


class PlotProperties(models.Model):
    chromosome_num = models.IntegerField(validators=[
        MaxValueValidator(22),
        MinValueValidator(1)
    ])
    percentage = models.DecimalField(max_digits=3, decimal_places=0, default=Decimal(0),
                                     validators=PERCENTAGE_VALIDATOR)
    min_index = models.IntegerField(validators=[
            MaxValueValidator(100000000),
            MinValueValidator(1)
        ])
    # max_index = models.IntegerField(validators=[
    #         MaxValueValidator(100000000),
    #         MinValueValidator(1)
    #     ])
