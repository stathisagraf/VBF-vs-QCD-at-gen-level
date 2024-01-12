import numpy as np

import matplotlib.pyplot as plt

import awkward as ak

import uproot

from collections import Counter

from hist_plotting import hist #

def y(pt, eta, m):
    return 1/2*np.log((np.sqrt(pt*pt*np.cosh(eta)*np.cosh(eta)+m*m)+pt*np.sinh(eta))/(np.sqrt(pt*pt*np.cosh(eta)*np.cosh(eta)+m*m)-pt*np.sinh(eta)))



def path_data(path,req_data):
    
    data = uproot.open(path)
    
    if req_data == "delta_eta":
        jet_eta = data['Events']['L1PuppiAK4Jets_eta'].array()
        delta_eta = []
        for dau in range(32):
            mask_eta = data['Events']['L1PuppiAK4Jets_dau{}_pt'.format(dau)].array()>-1
            
            delta_eta.append(ak.flatten(jet_eta[mask_eta]- data['Events']['L1PuppiAK4Jets_dau{}_eta'.format(dau)].array()[mask_eta]))
        data_set = delta_eta
        
    if req_data == "delta_phi":
        jet_phi = data['Events']['L1PuppiAK4Jets_phi'].array()
        delta_phi = []
        for dau in range(32):
            mask_phi = data['Events']['L1PuppiAK4Jets_dau{}_pt'.format(dau)].array()>-1
            delta_phi.append(ak.flatten(jet_phi[mask_phi]- data['Events']['L1PuppiAK4Jets_dau{}_phi'.format(dau)].array()[mask_phi]))
        data_set = delta_phi
    
    if req_data == "delta_y":
        
        jet_pt = data['Events']['L1PuppiAK4Jets_pt'].array()
        jet_eta = data['Events']['L1PuppiAK4Jets_eta'].array()
        jet_mass = data['Events']['L1PuppiAK4Jets_mass'].array()
        
        jet_y = y(jet_pt,jet_eta,jet_mass)
        
        delta_y = []
        for dau in range(32):
            mask_pt = data['Events']['L1PuppiAK4Jets_dau{}_pt'.format(dau)].array()>-1
            dau_pt = data['Events']['L1PuppiAK4Jets_dau{}_pt'.format(dau)].array()[mask_pt]
            dau_eta = data['Events']['L1PuppiAK4Jets_dau{}_eta'.format(dau)].array()[mask_pt]
            dau_mass = data['Events']['L1PuppiAK4Jets_dau{}_mass'.format(dau)].array()[mask_pt]
            delta_y.append(ak.flatten(jet_y[mask_pt]-y(dau_pt,dau_eta,dau_mass)))
        data_set = delta_y
    
    if req_data == "mass":
        dau_mass = []
        for dau in range(32):
            dau_mass.append(ak.flatten(data['Events']['L1PuppiAK4Jets_dau{}_mass'.format(dau)].array()))
        
        data_set = np.array(dau_mass)

    
    return data_set



vbf_path = 'perfNano_VBFHInv_PU200.l1ctlayer1.root'

mass = path_data(vbf_path,"mass").flatten()


plt.hist(mass, bins = 500)
plt.xlabel("Mass (GeV)")
plt.title('Puppi Candidates mass (AK4)')
plt.xlim(0,1.5)
plt.yscale('log')
plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Mass_plot_full_candidates_AK4.pdf')
plt.show()


delta_eta = ak.flatten(path_data(vbf_path,"delta_eta"))


fig,ax = plt.subplots()
alpha = hist(ax,delta_eta,1,51,False,'blue','Delta Eta',y_scale_log=False)
print(alpha[1][1]-alpha[1][0])

plt.title('Delta eta between jet axis and daughters (AK4)')
plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Delta_eta_plot_full_candidates_AK4.pdf')
plt.show()

fig,ax = plt.subplots()
alpha = hist(ax,delta_eta,1,51,False,'blue','Delta Eta',y_scale_log=True)
print(alpha[1][1]-alpha[1][0])

plt.title('Delta eta between jet axis and daughters (AK4)')
plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Delta_eta_plot_full_candidates_log_AK4.pdf')
plt.show()


delta_phi = ak.flatten(path_data(vbf_path,"delta_phi"))


fig,ax = plt.subplots()
beta = hist(ax,np.arcsin(np.sin(delta_phi)),1,51,False,'blue','Delta Phi',y_scale_log=False)
print(beta[1][1]-beta[1][0])
ax.set_xlim(-0.6,0.6)
plt.title('Delta phi between jet axis and daughters (AK4')
plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Delta_phi_plot_full_candidates_AK4.pdf')
plt.show()

fig,ax = plt.subplots()
beta = hist(ax,np.arcsin(np.sin(delta_phi)),1,51,False,'blue','Delta Phi',y_scale_log=True)
print(beta[1][1]-beta[1][0])

plt.title('Delta phi between jet axis and daughters (AK4)')
plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Delta_phi_plot_full_candidates_log_AK4.pdf')
plt.show()



delta_y = ak.flatten(path_data(vbf_path,"delta_y"))


fig,ax = plt.subplots()
beta = hist(ax,delta_y,1,51,False,'blue','Delta y',y_scale_log=False)
print(beta[1][1]-beta[1][0])

plt.title('Delta y between jet axis and daughters (AK4)')
plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Delta_y_plot_full_candidates_AK4.pdf')
plt.show()


fig,ax = plt.subplots()
beta = hist(ax,(delta_y),1,51,False,'blue','Delta y',y_scale_log=True)
print(beta[1][1]-beta[1][0])

plt.title('Delta y between jet axis and daughters (AK4)')
plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Delta_y_plot_full_candidates_log_AK4.pdf')
plt.show()