import numpy as np 
import pickle
import sys
import pyemma
from matplotlib import rc
import matplotlib as mpl
import matplotlib.pyplot as plt
import tol_colors as tc

hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

fig_wid = 10
fig_hig = 7

def rgb_to_hex(r, g, b):
    return ('#{:02X}{:02X}{:02X}').format(r, g, b)

def density_plot(prob_dens):
    min_bin = axis_lim[feature_name][0]
    max_bin = axis_lim[feature_name][1]
    fig,axs = plt.subplots(1,1,figsize=(fig_wid,fig_hig))
    bins = 50
    binsSD = np.linspace(min_bin,max_bin,51)
    averageSD = [(binsSD[i]+binsSD[i+1])/2 for i in range(bins)]

    plt.plot(averageSD,np.mean(prob_dens,axis=1)[0,:],label='HU-210',color=rgb_to_hex(int(255*0.88), int(255*0.47), int(255*0.90)),linewidth=3)
    plt.fill_between(averageSD, np.mean(prob_dens,axis=1)[0,:] + np.std(prob_dens,axis=1)[0,:],
                np.mean(prob_dens,axis=1)[0,:] - np.std(prob_dens,axis=1)[0,:], color=rgb_to_hex(int(255*0.88), int(255*0.47), int(255*0.90)), alpha=0.3)    #Error bar plot


    plt.plot(averageSD, np.mean(prob_dens,axis=1)[1,:], label='MDMB-FUBINACA', color=rgb_to_hex(int(255*0.02), int(255*0.38), int(255*0.67)),linewidth=3)
    plt.fill_between(averageSD, np.mean(prob_dens,axis=1)[1,:] + np.std(prob_dens,axis=1)[1,:],
                np.mean(prob_dens,axis=1)[1,:] - np.std(prob_dens,axis=1)[1,:], color=rgb_to_hex(int(255*0.02), int(255*0.38), int(255*0.67)), alpha=0.3)    #Error bar plot

    axs.set_xlim(axis_lim[feature_name])
    axs.set_ylim([0,0.8])

    axs.set_xticks(range(axis_lim[feature_name][0],axis_lim[feature_name][1]+1,2),range(axis_lim[feature_name][0],axis_lim[feature_name][1]+1,2),fontsize=22)
    axs.set_yticks(np.round(np.arange(0,0.9,0.2),1),np.round(np.arange(0,0.9,0.2),1),fontsize=22)
    
    plt.xlabel(axis_label[feature_name] + ' (\AA)', **hfont,fontsize=28)
    plt.ylabel('Probability Density', **hfont,fontsize=28)

    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.legend(fontsize=22)
    plt.tight_layout()
    plt.savefig('./Main_Figure_8B.png',transparent=False,dpi =300)

def density_cal(feature,weights):
    min_bin = axis_lim[feature_name][0]
    max_bin = axis_lim[feature_name][1]

    prob_density = np.empty(50)
    nSD, binsSD= np.histogram(feature, bins=50, range = (min_bin,max_bin),density=True, weights=weights)
    prob_density = nSD

    return prob_density

def weight_feature_calculation(method,ligand,bt):
    if method == 'TRAM':
        if ligand == 'nps':
            cluster = 700
        else:
            cluster = 800
        
        dtrajs = pickle.load(open('./tram_files/CB1_' + ligand + '_dtrajs_95_bt_' + bt + '.pkl','rb'))
        txx_dtrajs = np.concatenate(dtrajs)
        unique_clusters = np.sort(np.unique(txx_dtrajs))
        w = np.zeros(len(txx_dtrajs))
        tram_obj = pickle.load(open('./tram_files/CB1_' + ligand + '_95_bt_' + bt + '_obj.pkl','rb'))
        stat_dis = tram_obj.stationary_distribution_full_state
        del tram_obj

        for i in range(cluster):
            if i in unique_clusters:
                w[np.where(txx_dtrajs==i)[0]] = stat_dis[i]/len(np.where(txx_dtrajs==i)[0])
        
        return w

if __name__=='__main__': 
    feature_name = 'npxxy_rmsd'
    axis_label = {'npxxy_rmsd':'NPxxY RMSD'}
    axis_lim = {'npxxy_rmsd':[0,8]}

    feature_list = []
    weights_list = []
    prob_dens = np.empty([2,1,50])
    for i,ligand in enumerate(['nps','classical']):
        for j,bt in enumerate(['0']):
            feature = pickle.load(open('./feature_files/CB1_' + ligand + '_npxxy_rmsd_95_bt_' + bt + '.pkl','rb'))
            weights = weight_feature_calculation('TRAM',ligand,bt)
            prob_dens[i,j,:] = density_cal(feature,weights)
            del feature, weights 
    density_plot(prob_dens)

