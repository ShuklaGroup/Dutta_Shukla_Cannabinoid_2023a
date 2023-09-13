import os
import sys
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from matplotlib import rc
import seaborn as sns
import matplotlib as mpl

hfont = {'fontname':'Helvetica'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


def contact_freq_reader(state,bt,contact):
    resid = '324'
    filename = './frequency_files/CB1_nps_' + state + '_' + str(bt) + '_contact_freq.tsv'
    data = open(filename,'r').readlines()
    lig_str = 'X:LIG:' + resid
    for line in data:
        if line[0] == '#':
            continue
        lis = line.split()
        freq = float(lis[2])

        lis.remove(lig_str)
        temp = int(lis[0].split(':')[2])
        
        keys = list(contact[state].keys())

        if temp not in keys:
            contact[state][temp] = []
            contact[state][temp].append(freq)
        else:
            contact[state][temp].append(freq)

    return contact

def stable_interactions(contact,state):
    residues = list(contact[state].keys())
    stable_residues = []
    for resid in residues:
        if np.mean(contact[state][resid]) < 0.5:
            continue
        else:
            stable_residues.append(resid)
    return stable_residues

def heatmap(contact,stable_residues):
    df = pd.read_csv('./resid.csv',sep=',')
    new_contacts = {}
    residues_ticks = []
    for state in ['I2','Bound','I1','I3']:
        residues = list(contact[state].keys())
        new_contacts[state] = {}
        for ii,resid in enumerate(sorted(stable_residues)):
            if resid in residues:
                new_contacts[state][resid] = np.mean(contact[state][resid])
            else:
                new_contacts[state][resid] = 0
        
    for ii,resid in enumerate(sorted(stable_residues)):
        residues_ticks.append(str(df['resname'][int(resid)-1]) + str(df['OriNo'][int(resid)-1]) + '$^{' + str(df['Position'][int(resid)-1]) + '}$')

    df = pd.DataFrame.from_dict(new_contacts,orient='columns')
    fig, axs = plt.subplots(1,1,figsize=(5,12))
    ax = sns.heatmap(df,vmin=0,vmax=1.0,cbar=True,ax=axs,cmap=mpl.colormaps['Blues'],linewidths=1,square=True,linecolor='Black',cbar_kws={"shrink": 1.0,"fraction":0.15,"aspect":30,"label":'Contact Probability'})
    axs.tick_params(labelsize=30)
    axs.set_yticklabels(residues_ticks,rotation = 0,**hfont,fontsize=16)
    axs.set_xticklabels(['I2','Bound','I1','I3'],rotation = 30,**hfont,fontsize=16)
    plt.tight_layout()
    plt.savefig('./Main_Figure_6A.png',transparent=False,dpi =300)

if __name__=='__main__':
    stable_residues = []
    contact = {}
    for state in ['I2','Bound','I1','I3']:
        contact[state] = {}
        for bt in range(5):
            contact = contact_freq_reader(state,bt,contact)
        stable_residues.extend(stable_interactions(contact,state))
    
    heatmap(contact,list(set(stable_residues)))
