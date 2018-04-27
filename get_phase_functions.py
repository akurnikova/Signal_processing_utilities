# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 15:27:52 2016

@author: asya
"""

from __future__ import division
import numpy as np
from scipy.interpolate import interp1d

def get_phase_vals(vel_raw,L_filt,V_filt,fps,BrIn,resample_num):
    ## B is breathing array
    S = BrIn["valleys"]
    P = BrIn["peaks"]
    R = BrIn["rise"]
    L_sign = np.sign(np.diff(abs(L_filt)))
    V_sign = np.sign(np.diff(V_filt))
    
    for tp in np.arange(0,4):
        vel = np.copy(vel_raw)
        if tp == 0:
            vel[V_sign<0] = np.nan ## >0 = UP  <0 is DOWN
            vel[L_sign<0] = np.nan ## >0 is AWAY <0 is TO
        elif tp== 1:
            vel[V_sign<0] = np.nan ## >0 = UP  <0 is DOWN
            vel[L_sign>0] = np.nan ## >0 is AWAY <0 is TO    
        elif tp==2:
            vel[V_sign>0] = np.nan ## >0 = UP  <0 is DOWN
            vel[L_sign<0] = np.nan ## >0 is AWAY <0 is TO    
        elif tp==3:
            vel[V_sign>0] = np.nan ## >0 = UP  <0 is DOWN
            vel[L_sign>0] = np.nan ## >0 is AWAY <0 is TO    

        vel_phase_array = np.empty((len(R),resample_num*2))
        vel_phase_array.fill(np.nan)
        
        br_phase_array = np.empty((len(R),resample_num*2))
        br_phase_array.fill(np.nan)
        
        for i, i_start in enumerate(S):
            if (len(S)-2) < i:
                break
            i_end = S[i+1]
            i_mid = P[i+1]
            
            seg1 = np.linspace(-np.pi, 0, i_mid-i_start+1)
            seg2 = np.linspace(0,np.pi, i_end-i_mid)
            seg = np.concatenate((seg1[:-1],seg2))
            
            if len(seg) < 10:
                continue
            
            xp = np.linspace(-np.pi, np.pi, resample_num*2+1)    
            B = interp1d(seg, BrIn["breathe"][i_start:i_end])
            V = interp1d(seg, vel[i_start:i_end])
            
            vel_phase_array[i,:] = V(xp[:-1])
            br_phase_array[i,:] = B(xp[:-1])
        ## End iteration over breathe onsets
        
        #idx_notnan = ~np.isnan(vel_phase_array[:,0])
        idx_basal = (BrIn["instrate"][R] < 2)
        idx_sniff = (BrIn["instrate"][R] > 4)
    
           
        new_dict = {'vel_phase_array':vel_phase_array,
            'idx_sniff':idx_sniff,
            'idx_basal':idx_basal,
            'br_phase_array':br_phase_array,
            }
    
        if tp == 0:
            dat_UP_AWAY = new_dict
        elif tp == 1:
            dat_UP_TO = new_dict
        elif tp == 2:
            dat_DOWN_AWAY = new_dict
        elif tp == 3:
            dat_DOWN_TO = new_dict
        else:
            print "some sort of problem"
    ## end iteration over nose directions
    print "now return"
    return dat_UP_AWAY,dat_UP_TO,dat_DOWN_AWAY,dat_DOWN_TO
## end getter function

def get_phase_nosort(vel_raw,fps,BrIn,resample_num):
    ## B is breathing array
    S = BrIn["valleys"]
    P = BrIn["peaks"]
    R = BrIn["rise"]
    
    vel_phase_array = np.empty((len(R),resample_num*2))
    vel_phase_array.fill(np.nan)
        
    br_phase_array = np.empty((len(R),resample_num*2))
    br_phase_array.fill(np.nan)
        
    for i, i_start in enumerate(S):
        if (len(S)-2) < i:
            break
        i_end = S[i+1]
        i_mid = P[i+1]
        
        seg1 = np.linspace(-np.pi, 0, i_mid-i_start+1)
        seg2 = np.linspace(0,np.pi, i_end-i_mid)
        seg = np.concatenate((seg1[:-1],seg2))
        
        if len(seg) < 30:
            continue
        
        xp = np.linspace(-np.pi, np.pi, resample_num*2+1)    
        B = interp1d(seg, BrIn["breathe"][i_start:i_end],kind='linear')
        V = interp1d(seg, vel_raw[i_start:i_end],kind='linear')
        
        vel_phase_array[i,:] = V(xp[:-1])
        br_phase_array[i,:] = B(xp[:-1])
    ## End iteration over breathe onsets
        
    #idx_notnan = ~np.isnan(vel_phase_array[:,0])
    idx_basal = (BrIn["instrate"][R] < 2)
    idx_sniff = (BrIn["instrate"][R] > 4)
    
    
    new_dict = {'vel_phase_array':vel_phase_array,
        'idx_sniff':idx_sniff,
        'idx_basal':idx_basal,
        'br_phase_array':br_phase_array,
        }
    
    print "now return"
    return new_dict
## end getter function
    
def get_phase_fromonset(X_raw, V_filt, fps, BrIn, before_range,after_range):
    ## B is breathing array
    time_range = before_range+after_range
    R = BrIn["rise"]
    
    X_onset_array = np.empty((len(R),time_range*216))
    X_onset_array.fill(np.nan)
    
    br_onset_array = np.empty((len(R),time_range*216))
    br_onset_array.fill(np.nan)
    
    for i, r in enumerate(R):
        i_start = r - before_range*fps
        i_end = r + after_range*fps
        if i_start<0:
            continue
        if i_end>len(X_raw):
            break
    
    if fps == 1000:
        seg = np.linspace(before_range, after_range,time_range*fps)
        segN =np.linspace(before_range, after_range,time_range*216)
        B = interp1d(seg, BrIn["breathe"][i_start:i_end])
        X = interp1d(seg, X_raw[i_start:i_end])
        Bline = B(segN)
        Xline = X(segN)
    else:
        Bline = BrIn["breathe"][i_start:i_end]
        Xline = X_raw[i_start:i_end]
    
    #if len(Xline) < 10:
    #    continue
    
    X_onset_array[i,:] = Xline
    br_onset_array[i,:] = Bline
    ## End iteration over breathe onsets
    
    V_sign = np.sign(np.diff(V_filt))
    idx_basal = (BrIn["instrate"][R] < 2)
    idx_sniff = (BrIn["instrate"][R] > 4)
    idx_up = V_sign>0
    idx_dn = V_sign<0
    
    new_dict = {'X_onset_array':X_onset_array,
        'idx_sniff':idx_sniff,
        'idx_basal':idx_basal,
        'idx_up':idx_up,
        'idx_dn':idx_dn,
        'br_onset_array':br_onset_array,
        }

    return new_dict
## end getter function