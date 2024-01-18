import numpy as np
import matplotlib
import uproot
import awkward as ak
import matplotlib.pyplot as plt
import time
import sys
import vector
import numpy.ma as ma


def y(pt, eta, m):
    return 1/2*np.log((np.sqrt(pt*pt*np.cosh(eta)*np.cosh(eta)+m*m)+pt*np.sinh(eta))/(np.sqrt(pt*pt*np.cosh(eta)*np.cosh(eta)+m*m)-pt*np.sinh(eta)))



#the eta cut key differentiates the following cases: eta_cut_key==0 no eta cuts,  eta_cut_key==1 endcap region
# , eta_cut_key==-1 barrel region, eta_cut_key == 2  eta > leading_eta_cut, eta_cut_key == -2  eta < -leading_eta_cut
def pv(path, leading_pt_cut, subleading_pt_cut, leading_eta_cut, eta_cut_key, dau_pt_cuts = 0):
    
        
    data=uproot.open(path)
    
    
    number_jets_gen = data['Events']['nGenJets'].array()
    
    
    nonzero_gen = number_jets_gen > 1
    
    
    
    jet_pt_gen = data['Events']['GenJets_pt'].array()[nonzero_gen]
    
    
    sort_gen=ak.argsort(jet_pt_gen,ascending=False)
    leading_gen=sort_gen==0
    subleading_gen=sort_gen==1
    
    
   
    
    
    n_daus = 1
    
    
    n_daughters_l_gen = ak.flatten(data['Events']['GenJets_nDau'].array()[nonzero_gen][leading_gen] > n_daus)
    n_daughters_s_gen = ak.flatten(data['Events']['GenJets_nDau'].array()[nonzero_gen][subleading_gen] > n_daus)
    
    dau_mask_gen = n_daughters_l_gen*n_daughters_s_gen
   
    jet_pt_mask_leading = ak.flatten(data['Events']['GenJets_pt'].array()[nonzero_gen][leading_gen]>leading_pt_cut)
    jet_pt_mask_subleading = ak.flatten(data['Events']['GenJets_pt'].array()[nonzero_gen][subleading_gen]>subleading_pt_cut)
    
    
    if eta_cut_key==1:
        eta_cut_high = ak.flatten(data['Events']['GenJets_eta'].array()[nonzero_gen][leading_gen]>leading_eta_cut)
        eta_cut_low = ak.flatten(data['Events']['GenJets_eta'].array()[nonzero_gen][leading_gen]<-leading_eta_cut)
        eta_cut=eta_cut_high+eta_cut_low
        
        comb_gen_mask = dau_mask_gen*jet_pt_mask_leading*jet_pt_mask_subleading*eta_cut
        
    if eta_cut_key==-1:
        eta_cut_high = ak.flatten(data['Events']['GenJets_eta'].array()[nonzero_gen][leading_gen]<leading_eta_cut)
        eta_cut_low = ak.flatten(data['Events']['GenJets_eta'].array()[nonzero_gen][leading_gen]>-leading_eta_cut)
        eta_cut=eta_cut_high*eta_cut_low

        comb_gen_mask = dau_mask_gen*jet_pt_mask_leading*jet_pt_mask_subleading*eta_cut
        
    if eta_cut_key == 2:
        eta_cut_high = ak.flatten(data['Events']['GenJets_eta'].array()[nonzero_gen][leading_gen]>leading_eta_cut)
        comb_gen_mask = dau_mask_gen*jet_pt_mask_leading*jet_pt_mask_subleading*eta_cut_high
        
    if eta_cut_key == -2:
        eta_cut_low = ak.flatten(data['Events']['GenJets_eta'].array()[nonzero_gen][leading_gen]<-leading_eta_cut)
        comb_gen_mask = dau_mask_gen*jet_pt_mask_leading*jet_pt_mask_subleading*eta_cut_low


    if eta_cut_key == 0:
        comb_gen_mask = dau_mask_gen*jet_pt_mask_leading*jet_pt_mask_subleading
    
    n_dau_l_after_cuts_gen = ak.to_numpy(ak.flatten(data['Events']['GenJets_nDau'].array()[nonzero_gen][comb_gen_mask]))
    
    
    jet_pt_gen_leading = ak.to_numpy(ak.flatten(data['Events']['GenJets_pt'].array()[nonzero_gen][comb_gen_mask*leading_gen]))
    jet_pt_gen_subleading = ak.to_numpy(ak.flatten(data['Events']['GenJets_pt'].array()[nonzero_gen][comb_gen_mask*subleading_gen]))
    
    jet_eta_gen_leading = ak.to_numpy(ak.flatten(data['Events']['GenJets_eta'].array()[nonzero_gen][comb_gen_mask*leading_gen]))
    jet_eta_gen_subleading = ak.to_numpy(ak.flatten(data['Events']['GenJets_eta'].array()[nonzero_gen][comb_gen_mask*subleading_gen]))
    
    jet_phi_gen_leading = ak.to_numpy(ak.flatten(data['Events']['GenJets_phi'].array()[nonzero_gen][comb_gen_mask*leading_gen]))
    jet_phi_gen_subleading = ak.to_numpy(ak.flatten(data['Events']['GenJets_phi'].array()[nonzero_gen][comb_gen_mask*subleading_gen]))
    
    jet_mass_gen_leading = ak.to_numpy(ak.flatten(data['Events']['GenJets_mass'].array()[nonzero_gen][comb_gen_mask*leading_gen]))
    jet_mass_gen_subleading = ak.to_numpy(ak.flatten(data['Events']['GenJets_mass'].array()[nonzero_gen][comb_gen_mask*subleading_gen]))
    
    dau_pt_gen_leading = []
    dau_phi_gen_leading = []
    dau_eta_gen_leading =[]
    dau_mass_gen_leading = []
    
    dau_pt_gen_subleading = []
    dau_phi_gen_subleading = []
    dau_eta_gen_subleading =[]
    dau_mass_gen_subleading = []
    
    
    for dau in range(32):   
        dau_pt_gen_leading.append(ak.to_numpy(ak.flatten(data['Events']['GenJets_dau{}_pt'.format(dau)].array()[nonzero_gen][comb_gen_mask*leading_gen])))
        dau_eta_gen_leading.append(ak.to_numpy(ak.flatten(data['Events']['GenJets_dau{}_eta'.format(dau)].array()[nonzero_gen][comb_gen_mask*leading_gen])))
        dau_phi_gen_leading.append(ak.to_numpy(ak.flatten(data['Events']['GenJets_dau{}_phi'.format(dau)].array()[nonzero_gen][comb_gen_mask*leading_gen])))
        dau_mass_gen_leading.append(ak.to_numpy(ak.flatten(data['Events']['GenJets_dau{}_mass'.format(dau)].array()[nonzero_gen][comb_gen_mask*leading_gen])))
        
        dau_pt_gen_subleading.append(ak.to_numpy(ak.flatten(data['Events']['GenJets_dau{}_pt'.format(dau)].array()[nonzero_gen][comb_gen_mask*subleading_gen])))
        dau_eta_gen_subleading.append(ak.to_numpy(ak.flatten(data['Events']['GenJets_dau{}_eta'.format(dau)].array()[nonzero_gen][comb_gen_mask*subleading_gen])))
        dau_phi_gen_subleading.append(ak.to_numpy(ak.flatten(data['Events']['GenJets_dau{}_phi'.format(dau)].array()[nonzero_gen][comb_gen_mask*subleading_gen])))
        dau_mass_gen_subleading.append(ak.to_numpy(ak.flatten(data['Events']['GenJets_dau{}_mass'.format(dau)].array()[nonzero_gen][comb_gen_mask*subleading_gen])))
        
        
        
    y_jet_gen_leading = y(jet_pt_gen_leading, jet_eta_gen_leading, jet_mass_gen_leading)
    y_jet_gen_subleading = y(jet_pt_gen_subleading, jet_eta_gen_subleading, jet_mass_gen_subleading)
    
    sum_y_leading = 0
    sum_phi_leading = 0
    
    sum_y_subleading = 0 
    sum_phi_subleading = 0
    
    for dau in range(32):
        
        ntrue_dau_gen_leading=dau_pt_gen_leading[dau] < dau_pt_cuts
        
    
        dau_mask_pt_gen_leading=ma.masked_array(dau_pt_gen_leading[dau],mask=[ntrue_dau_gen_leading])
        dau_pt_complete_gen_leading=dau_mask_pt_gen_leading.filled(0)
        
        
        ntrue_dau_gen_subleading=dau_pt_gen_subleading[dau] < dau_pt_cuts
        
    
        dau_mask_pt_gen_subleading=ma.masked_array(dau_pt_gen_subleading[dau],mask=[ntrue_dau_gen_subleading])
        dau_pt_complete_gen_subleading=dau_mask_pt_gen_subleading.filled(0)
       
        if dau == 1:
            leading_dau_mask = ~ntrue_dau_gen_leading
            subleading_dau_mask = ~ntrue_dau_gen_subleading
        
        y_dau_gen_leading = y(dau_pt_complete_gen_leading, dau_eta_gen_leading[dau], dau_mass_gen_leading[dau])
        y_dau_gen_subleading = y(dau_pt_complete_gen_subleading, dau_eta_gen_subleading[dau], dau_mass_gen_subleading[dau])
        
        
        delta_y_leading = y_dau_gen_leading-y_jet_gen_leading
        delta_phi_leading = np.arcsin(np.sin(dau_phi_gen_leading[dau]-jet_phi_gen_leading))
        
        
        delta_y_subleading = y_dau_gen_subleading-y_jet_gen_subleading
        delta_phi_subleading = np.arcsin(np.sin(dau_phi_gen_subleading[dau]-jet_phi_gen_subleading))
        
        
        r_vec_leading = vector.array({
            "x" : delta_y_leading,
            "y" : delta_phi_leading})
        
        abs_r_leading = r_vec_leading.rho
        
        r_vec_subleading = vector.array({
            "x" : delta_y_subleading,
            "y" : delta_phi_subleading})
        
        abs_r_subleading = r_vec_subleading.rho
        
        
        kappa_vec_leading = dau_pt_complete_gen_leading*abs_r_leading
        
       
        
        sum_phi_leading += kappa_vec_leading*delta_phi_leading
        sum_y_leading += kappa_vec_leading*delta_y_leading
        
        
        kappa_vec_subleading = dau_pt_complete_gen_subleading*abs_r_subleading
        
        sum_y_subleading += kappa_vec_subleading*delta_y_subleading
        sum_phi_subleading += kappa_vec_subleading*delta_phi_subleading
        
        
        
    pull_vector_y_leading = sum_y_leading/jet_pt_gen_leading
    pull_vector_phi_leading = sum_phi_leading/jet_pt_gen_leading
    
    pull_vector_y_subleading = sum_y_subleading/jet_pt_gen_subleading
    pull_vector_phi_subleading = sum_phi_subleading/jet_pt_gen_subleading
    
    
    j_leading=vector.array({
                        "x" : y_jet_gen_leading[leading_dau_mask*subleading_dau_mask],
                        "y": jet_phi_gen_leading[leading_dau_mask*subleading_dau_mask]                
                        })
    
    j_subleading=vector.array({
                        "x" : y_jet_gen_subleading[leading_dau_mask*subleading_dau_mask],
                        "y": jet_phi_gen_subleading[leading_dau_mask*subleading_dau_mask]   
                        })
    
    pull_vector_leading = vector.array({
                        "x" : pull_vector_y_leading[leading_dau_mask*subleading_dau_mask],
                        "y": pull_vector_phi_leading[leading_dau_mask*subleading_dau_mask]
                        })
    
    pull_vector_subleading = vector.array({
                        "x" : pull_vector_y_subleading[leading_dau_mask*subleading_dau_mask],
                        "y": pull_vector_phi_subleading[leading_dau_mask*subleading_dau_mask]
                        })
    
    v12 = j_subleading-j_leading
    v21 = j_leading - j_subleading
    
    leading_relative_pa=v12.deltaphi(pull_vector_leading)
    subleading_relative_pa=v21.deltaphi(pull_vector_subleading)
    
    x=[leading_relative_pa,subleading_relative_pa,sum_phi_leading]   #put here every variable that you want to be outputted
        
    return pull_vector_y_leading[leading_dau_mask],pull_vector_phi_leading[leading_dau_mask],pull_vector_y_subleading[subleading_dau_mask],pull_vector_phi_subleading[subleading_dau_mask], x

# vbf_path = 'perfNano_VBFHInv_PU200.l1ctlayer1.root'
# qcd_path = 'perfNano_QCD_PU200.l1ctlayer1.root'


# vbf_pv_y_leading, vbf_pv_phi_leading, vbf_pv_y_subleading, vbf_pv_phi_subleading, x = pv(vbf_path,20,0,0,0)

# print(x[0])
# print(x[1])
# print(x[2])


# plt.hist([x[1],x[2]],bins = 100, histtype='step', color = ['black','red'])
# plt.yscale('log')
