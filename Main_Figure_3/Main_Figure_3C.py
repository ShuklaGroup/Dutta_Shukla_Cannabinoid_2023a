import matplotlib.pyplot as plt 
import glob
import pickle
import numpy as np
import tol_colors as tc
from matplotlib import rc
import matplotlib as mpl

hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

fig_wid = 10
fig_hig = 7

def rgb_to_hex(r, g, b):
    return ('#{:02X}{:02X}{:02X}').format(r, g, b)

def scatter_plot(Z_cor_per_traj,rmsd):
    fig,axs = plt.subplots(1,1,figsize=(fig_wid,fig_hig))
    print(rmsd[0].shape,rmsd[1].shape)
    plt.scatter(rmsd[0], Z_cor_per_traj[0], s=4, color='green', label='Replicate 1')
    plt.scatter(rmsd[1], Z_cor_per_traj[1], s=4, color='orange',label='Replicate 2')  
    axs.set_xlim([0,20])
    axs.set_ylim([0,30])

    axs.set_xlabel('RMSD (\AA)',fontsize=30)
    axs.set_ylabel('Z-component W$^{6.48}$ distance (\AA) ',fontsize=30)

    plt.xticks(np.arange(0,21,5),np.arange(0,21,5),fontsize=22)
    plt.yticks(np.arange(0,31,6),np.arange(0,31,6),fontsize=22)
    plt.legend(fontsize=18)
    plt.tight_layout()
    plt.savefig('./Main_Figure_3C.png',dpi=300)
    plt.close()



if __name__=='__main__':
    runs = [0,1] 
    Z_cor_per_traj = []
    rmsd = []    

    for i in runs:
        colvar_file = './colvar_files/CB1_classical_rep_' + str(i) +'_colvars.pkl'
        Z_cor_per_traj.append(pickle.load(open(colvar_file,'rb')))
        print(pickle.load(open(colvar_file,'rb')).shape)
        rmsd_file = './rmsd_files/CB1_classical_rep_' + str(i) + '_lig_rmsd.agr'
        f = np.loadtxt(rmsd_file,skiprows=8)
        rmsd.append(f[:,1]) 
        print(f.shape)

    scatter_plot(Z_cor_per_traj,rmsd)
