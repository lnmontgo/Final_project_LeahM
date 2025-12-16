# -*- coding: utf-8 -*-
"""
Created on Tue Dec  2 10:56:17 2025

@author: Leah Montgomery 

This script is the functions that I created for my ESCI 895 
final project. 

Functions: 
    1. read_csv_p: 
        this functions reads in the csv files for precipitation and returns a 
        dataframe with the appropriate precipitation data. 
    2. read_txt_dtw:
        this function reads in the txt file for the depth to water data. It 
        returns a dataframe with the depth to water for the aquifer. 
    3. correl_data:
        This function computes the cross correlation of precpitation data 
        against change in hydrolic head. 
    4. plot_corr:
        This function plots the correlations of the 4 aquifers that I am 
        using in this project. 
    
"""
#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate
import os
from pathlib import Path
from scipy.signal import savgol_filter

#%%
def read_csv_p(name, directory, subdirectory, start, end):
    """

    Parameters
    ----------
    name : string
        name of the csv file.
    directory : string
        The main directory name.
    subdirectory : string
        the subdirectory that the file name is located in.
    start : 
    end : 

    Returns
    -------
    data : DataFrame
        This is a dataframe with precipitation data.

    """
    data = pd.read_csv(directory/subdirectory/name, comment='#', header=0, parse_dates=['DATE'], index_col=['DATE'])
    data.replace(-99.9,np.nan, inplace=True)
    cols_list=data.columns.tolist()
    if 'PRCP' in cols_list:
        data['Precip'] = data['PRCP']
    data = data[['Precip']]
    data = data.resample('D').mean()
    data = data[pd.to_datetime(start):pd.to_datetime(end)]
    data = data.fillna(0)
    return data

#%%
def read_txt_dtw(name, directory, subdirectory, start, end):
    """

    Parameters
    ----------
    name : string
        name of the txt file.
    subdirectory : string
        the subdirectory that the file name is located in.

    Returns
    -------
    data : DataFrame
        This is a dataframe containing groundwater data.

    """
    data = pd.read_csv(directory/subdirectory/name, comment='#', header=[1], delimiter='\t',
                           parse_dates=['20d'])
    data.rename(columns={'20d':'Date'}, inplace=True)
    data.set_index('Date', inplace=True)
    cols_list = data.columns.tolist()
    if '14n' in cols_list and '14n.2' not in cols_list:
        data.rename(columns={'14n':'dtw'}, inplace=True)
    if '14n.2' in cols_list:
        data.rename(columns={'14n.2':"dtw"}, inplace=True)
    data = data[['dtw']]
    data = data.resample('D').mean()
    data = data[pd.to_datetime(start):pd.to_datetime(end)]
    data['diff'] = data['dtw'].diff()
    drop = data['diff'].abs() > 1.5
    data.loc[drop, 'dtw'] = np.nan
    data.interpolate(method='linear', limit_direction='both', inplace=True)
    
    return data

#%%
def correl_data(df_p, data_p, df_gw, data_gw):
    """

    Parameters
    ----------
    df_p : DataFrame
        Precipitation dataframe.
    df_dtw : DataFrame
        depth to water dataframe.
    data_p: string
        The name of the column with the desired Precip data 
    data_gw: string 
        the name of the column with the desired GW data 

    Returns
    -------
    corrs : DataFrame
        A dataframe with the correlations .
    lags : DataFrame
        A dataframe with the lag times.

    """
    p = df_p[data_p]
    gw = df_gw[data_gw]
    gw_detrended = gw - gw.rolling(60, center=True, min_periods=1).mean()
    lags = np.arange(0, 365)
    corrs = []
    for lag in lags:
        shifted = p.shift(lag).fillna(0)
        corr = shifted.corr(gw_detrended)
        corrs.append(corr)
    return corrs, lags

#%%
def plot_corr(p1, p2, p3, p4, gw1, gw2, gw3, gw4, data_p, data_gw):
    """

    Parameters
    ----------
    p1 : DataFrame
        First precipitation dataframe.
    p2 : DataFrame
        Second precipitation dataframe.
    p3 : DataFrame
        Third precipitation dataframe.
    p4 : DataFrame
        Fourth precipitation dataframe.
    gw1 : DataFrame
        First groundwater dataframe.
    gw2 : DataFrame
        Second groundwater dataframe.
    gw3 : DataFrame
        Third groundwater dataframe.
    gw4 : DataFrame
        Fourth groundwater dataframe.
    data_p : string
        Name of the column for precipitatio.
    data_gw : string
        Name of the column for groundwater data.

    Returns
    -------
    results : TYPE
        DESCRIPTION.

    """
    results = pd.DataFrame()
    
    in_corr, in_lags = correl_data(p1, data_p, gw1, data_gw)
    nh_corr, nh_lags = correl_data(p2, data_p, gw2, data_gw)
    co_corr, co_lags = correl_data(p3, data_p, gw3, data_gw)
    la_corr, la_lags = correl_data(p4, data_p, gw4, data_gw)
    
    results['IN'] = in_corr
    results['NH'] = nh_corr
    results['CO'] = co_corr
    results['LA'] = la_corr
    results['lags'] = in_lags
    
    fig2, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, sharex=True)
    fig2.set_size_inches(5,6)

    ax1.plot(nh_lags, nh_corr, color='blue', alpha=0.6)
    ax2.plot(co_lags, co_corr, color='red', alpha=0.6)
    ax3.plot(in_lags, in_corr, color='magenta', alpha=0.6)
    ax4.plot(la_lags, la_corr, color='green', alpha=0.6)
    ax1.set_title("New Hampshire")
    ax2.set_title("Colorado")
    ax3.set_title("Lafayette Indiana")
    ax4.set_title("Baton Rouge, Louisiana")
    ax3.set_ylabel('Correlation Coefficient')
    ax4.set_xlabel('Lags (days)')
    return results
    

#%%







