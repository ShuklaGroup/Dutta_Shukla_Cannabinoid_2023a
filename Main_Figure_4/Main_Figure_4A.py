import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rc
import pickle
import pyemma
import sys

hfont = {'fontname':'Helvetica'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def rgb_to_hex(r, g, b):
    return ('#{:02X}{:02X}{:02X}').format(r, g, b)

def bar_plot(Experiment,TRAM,MSM,TRAM_std,MSM_std):
    x1 = np.array([-0.1, 1.9, 3.9])
    x2 = np.array([0.1, 2.1, 4.1])

    fig, axs = plt.subplots(1,1,figsize=(10,7))
    
    axs.invert_yaxis()
    axs.bar(x1, MSM, yerr=MSM_std, color='Blue',width=0.2,label='MSM',alpha=0.8)
    axs.bar(x2, TRAM,yerr=TRAM_std,color='Orange', width=0.2,label='TRAM',alpha=0.8)

    axs.set_xlim([-1,5])
    x = np.arange(3.8,4.2,0.05)
    axs.plot(x,[Experiment[1]-Experiment[0]]*len(x),'--',c='Black',label='Experiment')
    axs.set_yticks(np.arange(-12,1,2),np.arange(-12,1,2), **hfont,fontsize=18)
    axs.set_yticks(np.arange(-12,1,2),np.arange(-12,1,2), **hfont,fontsize=18)

    axs.set_xticks([0,2,4],['HU-210', 'MDMB-Fubinaca', '$\Delta\Delta$G'], **hfont, fontsize=24)

    axs.set_ylabel('Binding Free Energy (kcal/mol)',**hfont,fontsize=28)

    plt.tight_layout()
    
    plt.legend(fontsize=18)

    plt.savefig('Main_Figure_4A.png',dpi=300)

    plt.close()

def free_energy_terms(txx,weights):
    R = 0.001987
    T = 300

    x_data = txx[:,0]*10
    y_data = txx[:,1]*10
    z_data = txx[:,2]*10

    x_data_min =  np.min(x_data)
    y_data_min =  np.min(y_data)
    x_data_max =  np.max(x_data)
    y_data_max =  np.max(y_data)
    z_data_min =  np.min(z_data)
    z_data_max =  np.max(z_data)

    x_hist_lim_low =  x_data_min -0.5
    y_hist_lim_low =  y_data_min -0.5
    x_hist_lim_high = x_data_max +0.5
    y_hist_lim_high = y_data_max  +0.5
    z_hist_lim_low =  z_data_min -0.5
    z_hist_lim_high = z_data_max  +0.5

    hist,edges = np.histogramdd(txx[:,:3]*10,bins=(x_bins,y_bins,z_bins),range = [[x_hist_lim_low,x_hist_lim_high],[y_hist_lim_low,y_hist_lim_high],[z_hist_lim_low,z_hist_lim_high]],density=True,weights=weights)

    xedge =  edges[0]
    yedge =  edges[1]
    zedge =  edges[2]

    x_bin_size = xedge[1]-xedge[0]
    y_bin_size = yedge[1]-yedge[0]
    z_bin_size = zedge[1]-zedge[0]

    free_energy = -R*T*np.log(hist*x_bin_size*y_bin_size*z_bin_size)
    min_free_energy= np.min(free_energy)
    delta_free_energy = free_energy - min_free_energy
    free_energy_flatten = delta_free_energy.flatten()
    del delta_free_energy,free_energy,hist
    max_free_energy = np.mean(sorted(free_energy_flatten[np.invert(free_energy_flatten==np.inf)])[-100:])
    binding_volume = np.sum(np.exp(free_energy_flatten[np.invert(free_energy_flatten==np.inf)]/(-R*T))*x_bin_size*y_bin_size*z_bin_size)
    volume_term = -R*T*np.log(binding_volume/1661)
    return volume_term - max_free_energy


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

    else:
        msm = pickle.load(open('./msm_files/CB1_' + ligand + '_bt_95_' + bt + '_msm.pkl','rb'))
        weights=np.concatenate(msm.trajectory_weights())

        return weights

if __name__=='__main__':
    R = 0.001987
    T = 300

    x_bins = 20
    y_bins = 20
    z_bins = 50

    Experiment = np.array([R*T*np.log(1.14*1e-9),R*T*np.log(0.061*1e-9)])
    Energy = np.empty([2,2,3])

    for i,method in enumerate(['TRAM','MSM']):
        for j,ligand in enumerate(['classical','nps']):
            for k,bt in enumerate(['0','1','2']):
                feature = pickle.load(open('./feature_files/CB1_classical_feature_' + method + '_95_bt_' + bt + '.pkl','rb'))
                weights = weight_feature_calculation(method,ligand,bt) 
                Energy[i,j,k] = free_energy_terms(feature,weights)        
                del feature,weights 

    TRAM_ = [np.mean(Energy[0,0,:]),np.mean(Energy[0,1,:])]
    TRAM_std = [np.std(Energy[0,0,:]),np.std(Energy[0,1,:])]

    MSM_ = [np.mean(Energy[1,0,:]),np.mean(Energy[1,1,:])]
    MSM_std = [np.std(Energy[1,0,:]),np.std(Energy[1,1,:])]

    ddG = []
    for i in range(3):
        for j in range(3):
            ddG.append(Energy[0,0,i] - Energy[0,1,j])

    TRAM_.append(np.mean(ddG))
    TRAM_std.append(np.std(ddG))

    ddG = []
    for i in range(3):
        for j in range(3):
            ddG.append(Energy[1,0,i] - Energy[1,1,j])

    MSM_.append(np.mean(ddG))
    MSM_std.append(np.std(ddG))

    bar_plot(Experiment,TRAM_,MSM_,TRAM_std,MSM_std)




