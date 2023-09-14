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
    bts = [25,50,75,95]
    nps_diss = np.empty([2,len(bts),3])
    clas_diss = np.empty([2,len(bts),3])

    for j,tech in enumerate(techs):
        macrostates = pickle.load(open('./cluster_index/CB1_nps_' + tech + '_macrostate_clusters.pkl','rb'))
        A = macrostates[0]
        B = macrostates[-1]
        for k,bt in enumerate(bts):
            nps_diss[j,k,:] = TPT_bt(A,B,bt,tech,'nps')

    for j,tech in enumerate(techs):
        macrostates = pickle.load(open('./cluster_index/CB1_classical_' + tech + '_macrostate_clusters.pkl','rb'))
        A = macrostates[1]
        B = macrostates[-1]
        for k,bt in enumerate(bts):
            clas_diss[j,k,:] = TPT_bt(A,B,bt,tech,'classical')

    nps_diss,clas_diss = nps_diss/10000,clas_diss/10000

    fig,axs = plt.subplots(1,1,figsize=(10,7))
    axs.set_yscale('symlog')
    count = 1
    for k,bt in enumerate(bts):
        diss_time = [[],[]]
        for j,tech in enumerate(techs):
            for ij in range(3):
                for ik in range(3):
                    diss_time[j].append(clas_diss[j,k,ij]-nps_diss[j,k,ik])    
        
        print(diss_time[0])
        box = plt.boxplot(diss_time,positions=[count,count+1],widths = 0.6)
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
    
    axs.set_xticks([1.5, 4.5, 7.5, 10.5], ['0.25','0.50','0.75','0.95'],**hfont, fontsize=24)
    axs.set_yticks([1e6,1e4,1e2,0,-1e2,-1e4,-1e6], ['10$^6$','10$^4$','10$^2$',0,'-10$^2$','-10$^4$','-10$^6$'],**hfont, fontsize=24)

    axs.legend([box["boxes"][0], box["boxes"][1]], ['MSM', 'TRAM'],loc='lower right',fontsize=24)
    plt.xlabel('Fraction of unbiased trajectories',**hfont, fontsize=30)
    plt.ylabel('$\Delta$ dissociation time  ($\mu$s)',**hfont, fontsize=30)
    
    plt.tight_layout()
    plt.savefig('Main_Figure_4D.png',dpi=300)
    plt.close()


