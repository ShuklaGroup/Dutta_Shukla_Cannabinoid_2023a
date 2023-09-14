import pyemma
import glob
import pickle
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np
import sys
from matplotlib import rc

hfont = {'fontname':'Helvetica'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def rgb_to_hex(r, g, b):
    return ('#{:02X}{:02X}{:02X}').format(r, g, b)

def update_index(A,active_set):
    new_index = []
    for i in A:
        if i in active_set:
            new_index.append(np.where(active_set==i)[0][0])

    return new_index

def TPT_bt(A,B,frac,tech,ligand):
    if tech == 'TRAM':
        filenames = glob.glob('./tram_files/CB1_' + ligand + '_' + str(frac) + '_bt_[0-2]_tram_msm_obj.pkl')
        timescale = []
        for filename in filenames:
            tram_msm =  pickle.load(open(filename,'rb'))
            active_set = tram_msm.active_set
            A_new = update_index(A,active_set)
            B_new = update_index(B,active_set)
            tpt = pyemma.msm.tpt(tram_msm,A_new,B_new)
            timescale.append(tpt.mfpt)
        return timescale

    else:
        filenames = glob.glob('./msm_files/CB1_' + ligand + '_bt_' + str(frac) + '_[0-2]_msm.pkl')
        timescale = []
        for filename in filenames:
            msm =  pickle.load(open(filename,'rb'))
            active_set = msm.active_set
            A_new = update_index(A,active_set)
            B_new = update_index(B,active_set)
            tpt = pyemma.msm.tpt(msm,A_new,B_new)
            timescale.append(tpt.mfpt)
        return timescale

if __name__=='__main__':
    techs = ['MSM','TRAM']
    frac = 95
    nps_diss = np.empty([2,3])
    nps_asso = np.empty([2,3])
    clas_diss = np.empty([2,3])
    clas_asso = np.empty([2,3])

    for j,tech in enumerate(techs):
        macrostates = pickle.load(open('./cluster_index/CB1_nps_' + tech + '_macrostate_clusters.pkl','rb'))
        A = macrostates[0]
        B = macrostates[-1]
        nps_diss[j,:] = TPT_bt(A,B,frac,tech,'nps')
        nps_asso[j,:] = TPT_bt(B,A,frac,tech,'nps')

    for j,tech in enumerate(techs):
        macrostates = pickle.load(open('./cluster_index/CB1_classical_' + tech + '_macrostate_clusters.pkl','rb'))
        A = macrostates[1]
        B = macrostates[-1]
        clas_diss[j,:] = TPT_bt(A,B,frac,tech,'classical')
        clas_asso[j,:] = TPT_bt(B,A,frac,tech,'classical')

    nps_diss,clas_diss = nps_diss/10000,clas_diss/10000
    nps_asso,clas_asso = nps_asso/10000,clas_asso/10000

    fig,axs = plt.subplots(1,1,figsize=(10,7))
    axs.set_yscale('log')
    count = 1
    box = plt.boxplot([clas_asso[0,:],clas_asso[1,:]],positions=[count,count+1],widths = 0.6)
    for item in ['whiskers', 'caps']:
        #print(item,box[item])
        plt.setp(box[item][0], color="Blue")
        plt.setp(box[item][1], color="Blue")
        plt.setp(box[item][2], color="Orange")
        plt.setp(box[item][3], color="Orange")

    for item in ['boxes', 'fliers','medians']:
        print(item,box[item])
        plt.setp(box[item][0], color="Blue")
        plt.setp(box[item][1], color="Orange")
    count += 3
    box = plt.boxplot([nps_asso[0,:],nps_asso[1,:]],positions=[count,count+1],widths = 0.6)
    print(nps_asso[0,:],nps_asso[1,:])
    for item in ['whiskers', 'caps']:
        #print(item,box[item])
        plt.setp(box[item][0], color="Blue")
        plt.setp(box[item][1], color="Blue")
        plt.setp(box[item][2], color="Orange")
        plt.setp(box[item][3], color="Orange")

    for item in ['boxes', 'fliers','medians']:
        print(item,box[item])
        plt.setp(box[item][0], color="Blue")
        plt.setp(box[item][1], color="Orange")

    axs.set_ylim([1e0,1e2])
    axs.set_xticks([1.5, 4.5], ['HU-210','MDMB-Fubinaca'],**hfont, fontsize=24)
    axs.set_yticks([1e2,1e1,1e0], ['10$^2$','10$^1$','10$^0$'],**hfont, fontsize=24)

    axs.legend([box["boxes"][0], box["boxes"][1]], ['MSM', 'TRAM'],loc='lower right',fontsize=24)
    plt.ylabel('Binding time ($\mu$s)',**hfont, fontsize=30)

    plt.tight_layout()
    plt.savefig('Main_Figure_4B.png',dpi=300)
    plt.close()



