import numpy as np 
import matplotlib.pyplot as plt 
import pyemma
import mdtraj as md
import os 
from sklearn.preprocessing import MinMaxScaler

def kl_divergence(p, q):
    kl_div = 0
    for i in range(len(p)):
        if not(p[i] == 0 or q[i] == 0):
            kl_div += p[i] * np.log2(p[i]/q[i])
    return kl_div


def kl_calculation(r,q):
    totdist = []
    ref_filename = './kl_div/CB1_nps_' + r + '_0_chv.npy'
    txx_ref = 1/np.load(ref_filename)

    query_filename = './kl_div/CB1_nps_' + q + '_0_chv.npy'
    txx_query = 1/np.load(query_filename)

    number_of_feature = txx_query.shape[1]

    kl_avg = np.empty([number_of_feature,1])

    for i in range(number_of_feature):
        x_data = txx_ref[:,i]
        y_data = txx_query[:,i]
        minimum = np.min(np.minimum(x_data,y_data))
        maximum = np.max(np.maximum(x_data,y_data))
        bins = np.arange(minimum,maximum,(maximum-minimum)/100)
        xhist,xedges = np.histogram(x_data,bins=bins,density=True)
        yhist,yedges = np.histogram(y_data,bins=bins,density=True)
        x_prob = xhist * np.diff(xedges)
        y_prob = yhist * np.diff(yedges)
        kl_avg[i,0] = (kl_divergence(x_prob,y_prob) + kl_divergence(y_prob,x_prob))/2

    np.save('./kl_div/CB1_nps_'+ r + '_' + q + '_kl_avg_score.npy',kl_avg)

    return kl_avg

def pdb_b_factor(r,q,j):
    f = open('./kl_div/CB1_nps_ref_' + q +'_0.pdb','r')
    fh = open('./kl_div/CB1_nps_ref_'+ r + '_' + q +'_b.pdb','w')

    data = f.readlines()
    b_values = {}
    for i in range(len(weights_per_residue)):
        temp = "{:.2f}".format(weights_per_residue[i,int(j)])
        b_values[str(indexs[i]+1)] = str(temp)

    print(b_values)
    for line in data:
        if line[:4] == 'ATOM' and not (int(line[23:26]) > 323):
            modified_line = line.replace('1.00  0.00','1.00  ' + b_values[str(int(line[23:26]))])
            fh.write(modified_line)
        else:
            fh.write(line)

    f.close()
    fh.close()

if __name__=='__main__':
    ref_state = 'I1'
    top = './kl_div/CB1_nps_strip.prmtop'
    a = np.load('./kl_div/CB1_nps_I1_0_chv_ind.npy')
    indexs = np.unique(a)
    weights_per_residue = np.empty([len(indexs),5])

    for i,state in enumerate(['I2','Bound','I1','I3']):
        filename = './kl_div/CB1_nps_' + state + '_0_chv.npy'
    
    for j,state2 in enumerate(['I2','Bound','I1','I3']):
        kl_file = './kl_div/CB1_nps_' + ref_state +  '_' + state2 + '_kl_avg_score.npy'
        if not os.path.exists(kl_file):
            kl_div = kl_calculation(ref_state,state2)
        
        kl_div = np.load(kl_file)


        for i in range(len(indexs)):
            position = np.where(a==indexs[i])[0]
            weights_per_residue[i,j] = np.sum(kl_div[position])
        
    
    min_max_scaler = MinMaxScaler()
    min_max_scaler.fit(weights_per_residue)
    weights_per_residue = min_max_scaler.transform(weights_per_residue)

    for j,state2 in enumerate(['I2','Bound','I1','I3']):
        pdb_b_factor(ref_state,state2,j)



    
