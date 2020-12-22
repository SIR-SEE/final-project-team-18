#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 13:08:03 2020

@author: annabeljerre
"""
import numpy as np
import matplotlib.pyplot as plt

#Variabler som kan varieras

N=10000 #populationen
β=2 #antal en infekterad person smittar varje dag
D=7 #antal dagar en persn är infekterad och kan smitta andra
γ = 1/D #andel av de infekterade som tilfrisknar per dag
R = β / γ #totalt antal smittade per person

#S(t): number of people susceptible on day t
#I(t): number of people infected on day t
#R(t): number of people recovered on day t

