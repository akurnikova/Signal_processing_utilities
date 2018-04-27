# -*- coding: utf-8 -*-
"""
Processing for breathing signal for mouse (implanted thermistor)
Stacy 2017
"""
from __future__ import division

import scipy.signal as sig
from utils_filtering import util_lowpass, util_hipass
import numpy as np


def util_proc_breathe(dat, fps):
    br_sig = util_get_breathe(dat,fps)
    pk, val = util_get_breathe_peaks(br_sig,fps)
    rise = util_calc_breathe_pct(br_sig, pk, val, 10)
    instrate, duration = util_get_br_rate(br_sig, rise, fps)
    breathe = {'breathe':br_sig, 'peaks':pk,'valleys':val,'rise':rise,'duration':duration,'instrate':instrate}
    return breathe


def util_get_breathe(dat,fps):
    ts = np.shape(dat['breathe_raw'])
    if len(ts)==1:
        br = dat['breathe_raw']
    else:
        if ts[1]>ts[0]:
            br = dat['breathe_raw'][0,:]
        else:
            br = dat['breathe_raw'][:,0]
    br_lp = util_lowpass(br-np.mean(br),fps,30,3)
    br_hp = util_hipass(br_lp,fps,1,1)
    br_sig = br_hp-util_lowpass(br_hp,fps,0.1,1)
    return br_sig
    
    
def util_get_breathe_peaks(br_lp,fps):
    print 'starting hilbert'
    downsample = 1
    if fps >=10000:
        downsample = 10
    br_lp = br_lp[::downsample] 
    br_ang = np.angle(sig.hilbert(-br_lp)) ## this step is slow
    br_reset = np.diff(br_ang)
    br_pk = np.where(br_reset < -np.pi)[0]
    print 'end hilbert'
    return downsample*br_pk, downsample*calc_breathe_vals(br_lp, br_pk)
  
    
def calc_breathe_vals(br_lp, pk):
    # Recalculate peaks as maxima between valleys 
    br_val_new = np.zeros(len(pk)-1) # Initialize
    for i, p in enumerate(pk[:-1]):
        br_seg1 = -br_lp[p:pk[i+1]]
        signal_peak = np.amax(br_seg1)
        br_peak_time2 = (p + np.where(br_seg1 == signal_peak))[0,0]
        br_val_new[i] = br_peak_time2
    return br_val_new.astype(int)

    
def util_calc_breathe_pct(br_sig, pk, val, pct):
    pk = pk[1:]
    max_btw_pks = br_sig[pk]-br_sig[val]
    thresh_val = br_sig[val]+(pct/100)*(max_btw_pks)  
    risetimes_xpct = np.zeros(len(thresh_val))

    for ievent,t in enumerate(thresh_val):
        rise = br_sig[val[ievent]:pk[ievent]]
        if np.any(rise>=t):
            thresh = np.where(rise>=t)[0][0]
        else:
            thresh = 0
        risetimes_xpct[ievent] = val[ievent]+thresh;
    return risetimes_xpct.astype(int)


def util_get_br_rate(br_sig, rise, fps):
    instrate = np.zeros(len(br_sig))
    duration = np.diff(rise)/fps
    instrate[0:rise[0]] = 1/duration[1]

    for i, dur in enumerate(duration):
        i_start = rise[i]
        i_end = rise[i+1]
        instrate[i_start:i_end] = 1/dur;
    return instrate, duration