# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 19:51:19 2016

@author: asya
"""
from __future__ import division
import numpy as np
import scipy.signal as sig
from scipy import stats as st
import matplotlib.pyplot as plt

def util_lowpass(signal_raw,fps,hicut,order):
    mean_calc = np.mean(signal_raw) # subtract and add back mean
    S = signal_raw - mean_calc
    b, a = sig.butter(order,hicut/fps,'lowpass')
    signal_filt = sig.filtfilt(b,a,S)
    return signal_filt + mean_calc
    
def util_hipass(signal_raw,fps,cut,order):
    mean_calc = np.mean(signal_raw) # subtract and add back mean
    S = signal_raw - mean_calc
    b, a = sig.butter(order,cut/fps,'highpass')
    signal_filt = sig.filtfilt(b,a,S)
    return signal_filt + mean_calc
    
def get_angular_mean(r,ang,confidence=0.95,plotvecs = 1,limplt = 0):
    ## RUN AS mX,cuX,clX, mY,cuY,clY = get_angular_mean(r,ang)
    n = len(r)
    X_coords = r*np.cos(ang)
    Y_coords = r*np.sin(ang)
    mX,seX = np.nanmean(X_coords), st.sem(Y_coords)
    mY,seY = np.nanmean(Y_coords), st.sem(Y_coords)
    hX = seX * st.t._ppf((1+confidence)/2., n-1)
    hY = seY * st.t._ppf((1+confidence)/2., n-1)
    mX,cuX,clX, mY,cuY,clY = mX, mX + hX,mX - hX, mY,mY + hY,mY - hY
    if plotvecs == 1:
        plt.plot([np.arctan2(mY,mX),np.arctan2(mY,mX)],[limplt,(mY**2+mX**2)**0.5],'g')
        plt.plot([np.arctan2(clY,clX),np.arctan2(cuY,cuX)],[(cuX**2+cuY**2)**.5,(cuX**2+cuY**2)**.5],'g')
        plt.plot([np.arctan2(cuY,cuX),np.arctan2(clY,clX)],[(clX**2+clY**2)**.5,(clX**2+clY**2)**.5],'g')
        plt.plot([np.arctan2(clY,clX),np.arctan2(clY,clX)],[(cuX**2+cuY**2)**.5,(clX**2+clY**2)**.5],'g')
        plt.plot([np.arctan2(cuY,cuX),np.arctan2(cuY,cuX)],[(clX**2+clY**2)**.5,(cuX**2+cuY**2)**.5],'g') 
    return mX, mX + hX,mX - hX, mY,mY + hY,mY - hY
    
def get_angular_mean_nonzero(r,ang,thresh = 0.01,confidence=0.95,plotvecs = 1,limplt = 0):
    ## RUN AS mX,cuX,clX, mY,cuY,clY = get_angular_mean(r,ang)
    n = len(r)
    idk_keep = r > thresh
    X_coords = r[idk_keep]*np.cos(ang[idk_keep])
    Y_coords = r[idk_keep]*np.sin(ang[idk_keep])
    mX,seX = np.nanmean(X_coords), st.sem(Y_coords)
    mY,seY = np.nanmean(Y_coords), st.sem(Y_coords)
    hX = seX * st.t._ppf((1+confidence)/2., n-1)
    hY = seY * st.t._ppf((1+confidence)/2., n-1)
    mX,cuX,clX, mY,cuY,clY = mX, mX + hX,mX - hX, mY,mY + hY,mY - hY
    if plotvecs == 1:
        plt.plot([np.arctan2(mY,mX),np.arctan2(mY,mX)],[limplt,(mY**2+mX**2)**0.5],'g',linewidth = 2)
        plt.plot([np.arctan2(clY,clX),np.arctan2(cuY,cuX)],[(cuX**2+cuY**2)**.5,(cuX**2+cuY**2)**.5],'g')
        plt.plot([np.arctan2(cuY,cuX),np.arctan2(clY,clX)],[(clX**2+clY**2)**.5,(clX**2+clY**2)**.5],'g')
        plt.plot([np.arctan2(clY,clX),np.arctan2(clY,clX)],[(cuX**2+cuY**2)**.5,(clX**2+clY**2)**.5],'g')
        plt.plot([np.arctan2(cuY,cuX),np.arctan2(cuY,cuX)],[(clX**2+clY**2)**.5,(cuX**2+cuY**2)**.5],'g') 
    return mX, mX + hX,mX - hX, mY,mY + hY,mY - hY
    
    
def average_by_angular_coords_nonzero(r,ang,thresh = 0.01,confidence=0.95,plotvecs = 1,limplt = 0, bins = 30):
    ## RUN AS mX,cuX,clX, mY,cuY,clY = get_angular_mean(r,ang)
    ## First exclude values undder threshold r
    idk_keep = r > thresh
    R = r[idk_keep]
    A = ang[idk_keep]
    
    R = R[~np.isnan(A)]
    A = A[~np.isnan(A)]
    ## Create average angle array
    A_av = np.linspace(np.min(A),np.max(A),bins+1)
    A_av = A_av[:-1]
    ## Initialize R averaging arrays
    binstep = 2*np.pi/bins
    R_av = np.zeros((bins,1))
    R_sem = np.zeros((bins,1))
    
    ## Loop over the angles
    for i, a1 in enumerate(A_av):
        idx_next = np.where((A>a1)*(A<=a1+binstep))[0]
        R_av[i] = np.mean(R[idx_next])
        R_sem[i] = np.std(R[idx_next])/(len(idx_next))**0.5
        
    R_sem = R_sem[:,0]
    R_av = R_av[:,0]

    if plotvecs == 1:
        plt.plot(A_av,R_av,'r',linewidth = 2)
        plt.fill_between(A_av,R_av+R_sem,R_av-R_sem,alpha = 0.3)
    return A_av, R_av, R_sem