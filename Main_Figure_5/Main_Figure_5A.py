import numpy as np
import matplotlib.pyplot as plt
import pickle
import pyemma
from matplotlib import rc
import matplotlib as mpl
import matplotlib.pyplot as plt
from multiprocessing import Pool
import tol_colors as tc

hfont = {'fontname':'Helvetica','fontweight':'bold'}
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

fig_wid = 10
fig_hig = 7

def two_D_free_energy(x_data,y_data,weights):
    x_bins = 200
    y_bins = 200
    R = 0.001987
    T = 300

    fig_wid = 10
    fig_hig = 7
    cmap = tc.tol_cmap('rainbow_PuBr')

    Max_energy = 10
    x_data_min =  np.min(x_data)
    y_data_min =  np.min(y_data)
    x_data_max =  np.max(x_data)
    y_data_max =  np.max(y_data)

    x_hist_lim_low =  x_data_min -0.5
    y_hist_lim_low =  y_data_min -0.5
    x_hist_lim_high = x_data_max +0.5
    y_hist_lim_high = y_data_max  +0.5

    x_lim_low = (int(np.min(x_data)/5))*5.0
    y_lim_low = (int(np.min(y_data)/5))*5.0
    x_lim_high = (int(np.max(x_data)/5) + 1)*5.0
    y_lim_high = (int(np.max(y_data)/5) + 1)*5.0

    hist= np.histogram2d(x_data,y_data, bins=[x_bins,y_bins],
                         range = [[x_hist_lim_low,x_hist_lim_high],[y_hist_lim_low,y_hist_lim_high]],
                         density= True,weights=weights)

    prob_density = hist[0]
    xedge = hist[1]
    yedge = hist[2]

    x_bin_size = xedge[1]-xedge[0]
    y_bin_size = yedge[1]-yedge[0]

    free_energy = -R*T*np.log(prob_density*x_bin_size*y_bin_size)
    min_free_energy= np.min(free_energy)

    delta_free_energy = free_energy - min_free_energy

    xx = [(xedge[i]+xedge[i+1])/2 for i in range(len(xedge)-1)]
    yy = [(yedge[i]+yedge[i+1])/2 for i in range(len(yedge)-1)]

    X, Y = np.meshgrid(xx,yy)
    fig, axs = plt.subplots(1,1,figsize=(fig_wid,fig_hig))

    cd =axs.contourf(xx,yy,delta_free_energy.T, np.linspace(0,Max_energy,Max_energy*5+1),
                     vmin=0.0, vmax=Max_energy,cmap=cmap)
    cbar = fig.colorbar(cd,ticks=range(Max_energy+1))
    cbar.ax.set_yticklabels(range(Max_energy+1),fontsize=22)                                        #ticklabels of color bar 
    cbar.ax.set_ylabel('Free Energy (Kcal/mol)', labelpad=15,**hfont,fontsize=24)

    plt.contour(X,Y,delta_free_energy.T, linewidths=0.25, levels=range(0,11,2), colors='black')
    axs.set_xlim([0,30])
    axs.set_ylim([0,30])

    axs.set_xticks(range(0,30+1,5),range(0,30+1,5),**hfont,fontsize=22)
    axs.set_yticks(range(0,30+1,5),range(0,30+1,5),**hfont,fontsize=22)
    

    plt.xlabel('TM5(W279$^{5.43}$) and Ligand distance (\AA)', **hfont,fontsize=26)
    plt.ylabel('TM7(S383$^{7.39}$) and Ligand distance (\AA)', **hfont,fontsize=26)

    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.tight_layout()

    plt.savefig('./Main_Figure_5A.png',transparent=False,dpi =300)
    

if __name__=='__main__':       
    data = pickle.load(open('Main_Figure_5A.pkl','rb'))

    x_data = data[0]
    y_data = data[1]

    dtrajs = pickle.load(open('CB1_nps_dtrajs.pkl','rb'))
    txx_dtrajs = np.concatenate(dtrajs)
    unique_clusters = np.unique(txx_dtrajs)
    number_per_uni_clus = [len(np.where(txx_dtrajs==i)[0]) for i in unique_clusters]

    tram_obj = pickle.load(open('CB1_nps_tram_obj.pkl','rb'))
    stat_dis = tram_obj.stationary_distribution_full_state

    weights =  np.array([stat_dis[txx_dtrajs[i]]/number_per_uni_clus[txx_dtrajs[i]] for i in range(len(txx_dtrajs))])
    two_D_free_energy(x_data,y_data,weights)
