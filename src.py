# -*- coding: utf-8 -*-
import os

import numpy as np
import pandas as pd
import scipy.ndimage
from matplotlib import pyplot as plt
from scipy import optimize

base_dir = os.path.dirname(__file__)
    
def experimental_data_read(filename):
    file_path = os.path.join(base_dir, filename)
    df = pd.read_excel(file_path)
    return df
    
def linear_region(df, y_start, y_stop):
    x = df["X"]
    y = df["Y"]
    
    df_linear = df[(y > y_start) & (y < y_stop)].iloc[:1000]
    return df_linear

def residuals(p, y, x):
    k, b = p
    return y - (k*x + b)
        
def curve_plot(i, df, df_linear):
    x = df["X"]
    y = df["Y"]
    x_linear = df_linear["X"]
    y_linear = df_linear["Y"]
  
    plt.figure()
    plt.plot(x, y, label = "experimental data")
    plt.plot(x_linear, y_linear, label = "linear part")
    plt.legend()
    plt.savefig(str(i), format="png")
    
def least_square_fitting(rediduals, p0, df_linear):
    x_linear = df_linear["X"]
    y_linear = df_linear["Y"]
    
    r = optimize.leastsq(residuals, p0, args=(y_linear, x_linear))
    k, b = r[0]
    return k

y_start  = float(raw_input("Enter y_start:"))
y_stop = float(raw_input("Enter y_stop:"))
k = []

for i in range(1,6):
    filename = '.'.join((str(i),'xls'))
    df = experimental_data_read(filename)
    df_linear = linear_region(df, y_start, y_stop)  
     
    curve_plot(i, df, df_linear)
    
    p0 = [1, 0]

    k.append(least_square_fitting(residuals, p0, df_linear))
    
slope = pd.Series(k)
slope.to_csv('slope.csv')
    

    
