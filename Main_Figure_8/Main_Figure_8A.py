import numpy as np
import os
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

def density_cal(feature,weights):
    min_bin = axis_lim[feature_name][0]
    max_bin = axis_lim[feature_name][1]

    prob_density = np.empty(bins)
    nSD, binsSD= np.histogram(feature, bins=bins, range = (min_bin,max_bin),density=True, weights=weights)
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
        del tram_obj,dtrajs

        for i in range(cluster):
            if i in unique_clusters:
                w[np.where(txx_dtrajs==i)[0]] = stat_dis[i]/len(np.where(txx_dtrajs==i)[0])

        return w

def bar_plot(prob_dens):
    fig,axs = plt.subplots(1,1,figsize=(10,7))
    random_number = np.random.normal(0,0.01,3)

    axs.bar([1,2],np.mean(prob_dens,axis=1),yerr=np.std(prob_dens,axis=1),width=0.4,label=['MDMB-FUBINACA','HU-210'],color=[rgb_to_hex(int(255*0.02), int(255*0.38), int(255*0.67)),rgb_to_hex(int(255*0.88), int(255*0.47), int(255*0.90))])
    axs.scatter(random_number+1,prob_dens[0,:],s=4,c='Black')
    axs.scatter(random_number+2,prob_dens[1,:],s=4,c='Black')
    axs.set_xlim([0.5,2.5])
    axs.set_ylim([0,0.15])
    #axs.set_ylabel('T$^{3.46}$-Y$^{5.58}$-Y$^{7.53}$ triad',**hfont,fontsize=28)
    axs.set_xticks([1,2],['MDMB-FUBINACA','HU-210'],**hfont, fontsize=20)
    axs.set_yticks([0,0.05,0.10,0.15],[0,0.05,0.10,0.15],**hfont, fontsize=24)
    axs.set_ylabel('Triad (T$^{3.46}$-Y$^{5.58}$-Y$^{7.53}$) Probability',**hfont,fontsize=24)

    plt.savefig('./Main_Figure_8A.png',dpi=300)
    plt.tight_layout()
    plt.close()
        


if __name__=='__main__': 
    bins = 50
    feature_name = 'TM35'
    ind = {'TM36':0,'TM57':1,'TG':3,'TM35':4,'TM37':5}
    axis_label = {'TM36':'TM3-TM6 distance','TM57':'TM5-TM7 distance','TM35':'TM3-TM5 distance','TM37':'TM3-TM7 distance'}
    axis_lim = {'TM36':[8,20],'TM57':[2,20],'TM35':[2,20],'TM37':[2,20]}

    feature_list = []
    weights_list = []
    prob_dens = np.empty([2,3])
    for i,ligand in enumerate(['nps','classical']):
        for j,bt in enumerate(['0','1','2']):
            temp = 0
            for k,feature_name in enumerate(['TM57','TM35','TM37']):
                feature = pickle.load(open('./feature_files/CB1_' + ligand + '_' + feature_name + '_95_bt_' + bt + '.pkl','rb')) 
                temp += np.where(feature<5,1,0)

            weights = weight_feature_calculation('TRAM',ligand,bt) 
            prob_dens[i,j] = np.sum(weights[temp==3])
            #print(ligand,bt,np.sum(np.multiply(feature,weights)),np.sum(weights[feature<6]))
            #del feature, weights
    bar_plot(prob_dens)
    #density_plot(prob_dens)

