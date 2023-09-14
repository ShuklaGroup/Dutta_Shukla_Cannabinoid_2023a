import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import pickle
import sys
import pandas as pd
from matplotlib import rc

hfont = {'fontname':'Helvetica'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def rgb_to_hex(r, g, b):
    return ('#{:02X}{:02X}{:02X}').format(r, g, b)


def bar_plot(nps_weights,classical_weights):
    df = pd.read_csv('./resid.csv',sep=',')
    fig,axs = plt.subplots(1,1,figsize=(12,7))
    random_number = np.random.normal(0,0.01,3)
    
    width = 0.4
    x1 = np.array([i - 0.2 for i in range(len(binding_pocket_residue))])
    x2 = np.array([i + 0.2 for i in range(len(binding_pocket_residue))])
    
    axs.bar(x1,np.mean(nps_weights,axis=1),yerr=np.std(nps_weights,axis=1),width=0.4,label='MDMB-FUBINACA',color=rgb_to_hex(int(255*0.02), int(255*0.38), int(255*0.67)))
    axs.bar(x2,np.mean(classical_weights,axis=1),yerr=np.std(classical_weights,axis=1),width=0.4,label='HU-210',color=rgb_to_hex(int(255*0.88), int(255*0.47), int(255*0.90)))
    
    temp = x1.reshape((len(x1),1)) + random_number
    axs.scatter(temp.reshape(-1),nps_weights.reshape(-1),s=4,c='Black')

    temp = x2.reshape((len(x2),1)) + random_number
    axs.scatter(temp.reshape(-1),classical_weights.reshape(-1),s=4,c='Black')

    axs.set_xlim([-1,len(binding_pocket_residue)])
    axs.set_ylim([0,0.3])
    
    x = np.arange(0,len(binding_pocket_residue),1)
    residues_ticks = []
    for ii,resid in enumerate(sorted(binding_pocket_residue)):
        residues_ticks.append(str(df['resname'][int(resid)]) + str(df['OriNo'][int(resid)]) + '$^{' + str(df['Position'][int(resid)]) + '}$')

    axs.set_xticks(x,residues_ticks,**hfont, fontsize=24,rotation=30)
    axs.set_yticks(np.round(np.arange(0,0.31,0.05),2),np.round(np.arange(0,0.31,0.05),2), **hfont,fontsize=18)

    axs.set_ylabel('Allosteric Weight',**hfont,fontsize=28)
    plt.legend(fontsize=18)
    plt.tight_layout()
    plt.savefig('./Main_Figure_8C.png',dpi=300)
    plt.close()
    


if __name__=='__main__':
    binding_pocket_residue = [81, 84, 85, 88, 104, 108, 187, 190, 294]

    NPXXY = [i for i in range(304,309)]

    count = 0
     
    nps_weights = np.empty([len(binding_pocket_residue),3])
    classical_weights = np.empty([len(binding_pocket_residue),3])

    for i in range(3):
        nps = pickle.load(open('./weights/CB1_nps_' + str(i) + '_edges_results_visual.pkl','rb'))
        classical = pickle.load(open('./weights/CB1_classical_' + str(i) + '_edges_results_visual.pkl','rb'))
        for j,bpi in enumerate(binding_pocket_residue):
            temp_nps = 0
            temp_classical = 0
            for npx in NPXXY:
                temp_nps += nps[bpi,npx]
                temp_classical += classical[bpi,npx]
        
            nps_weights[j,i] = temp_nps/len(NPXXY)
            classical_weights[j,i] = temp_classical/len(NPXXY)
    
    print(np.mean(nps_weights,axis=1),np.std(nps_weights,axis=1))
    print(np.mean(classical_weights,axis=1),np.std(classical_weights,axis=1))
    bar_plot(nps_weights,classical_weights)

