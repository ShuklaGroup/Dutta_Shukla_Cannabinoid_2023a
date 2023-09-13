import pyemma
import glob
import pickle
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np
import sys
from matplotlib import rc

def update_index(A,active_set):
    new_index = []
    for i in A:
        if i in active_set:
            new_index.append(np.where(active_set==i)[0][0])

    return new_index

def TPT_bt(A,B):
    filenames = glob.glob('./tram_files/CB1_classical_95_bt_[0-2]_tram_msm_obj.pkl')
    timescale = [] 
    for filename in filenames:
        tram_msm =  pickle.load(open(filename,'rb'))
        active_set = tram_msm.active_set
        A_new = update_index(A,active_set)
        B_new = update_index(B,active_set)
        tpt = pyemma.msm.tpt(tram_msm,A_new,B_new)
        timescale.append(tpt.mfpt)
    return np.mean(timescale),np.std(timescale)


if __name__=='__main__':
    timescales = []
    macro_ind = {'I1':0,'Bound':1,'I2':2,'I3':3,'Unbound':-1}
    Intial_state = ['Bound','I1','Bound','I2','I2','I3','I2','I3','Unbound','Unbound']
    Final_state = ['I1','Bound','I2','Bound','I3','I2','Unbound','Unbound','I2','I3']

    for I,F in zip(Intial_state,Final_state):
        macrostates = pickle.load(open('./CB1_classical_TRAM_macrostate_clusters.pkl','rb'))
        A = macrostates[macro_ind[I]]
        B = macrostates[macro_ind[F]]

        temp = np.empty(2)

        mean,std = TPT_bt(A,B)
        temp[0] = mean/10000
        temp[1] = std/10000
    
        print(I + ' ' + F + ' ' + str(temp[0]) + '+/-' + str(temp[1]))
