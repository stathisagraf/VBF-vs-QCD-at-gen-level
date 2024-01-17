import numpy as np
import matplotlib
import uproot
import awkward as ak
import matplotlib.pyplot as plt
import time
import sys
import vector


def masks_Lj():
    number_jets_gen = vbf_data['Events']['nGenJets'].array()
    number_jets_AK4 = vbf_data['Events']['nL1PuppiAK4Jets'].array()

    nonzero_gen = number_jets_gen > 1
    nonzero_AK4 = number_jets_AK4 > 1


    jet_pt_gen = vbf_data['Events']['GenJets_pt'].array()[nonzero_gen]
    jet_pt_AK4 = vbf_data['Events']['L1PuppiAK4Jets_pt'].array()[nonzero_AK4]

    sort_gen=ak.argsort(jet_pt_gen,ascending=False)
    leading_gen=sort_gen==0
    subleading_gen=sort_gen==1


    sort_AK4=ak.argsort(jet_pt_AK4,ascending=False)
    leading_AK4=sort_AK4==0
    subleading_AK4=sort_AK4==1

    n_daus = 1

    n_daughters_l_gen = vbf_data['Events']['GenJets_nDau'].array()[nonzero_gen][leading_gen] > n_daus
    n_daughters_s_gen = vbf_data['Events']['GenJets_nDau'].array()[nonzero_gen][subleading_gen] > -1



    dau_mask_gen = n_daughters_l_gen*n_daughters_s_gen



    n_daughters_l_AK4 = vbf_data['Events']['L1PuppiAK4Jets_nDau'].array()[nonzero_AK4][leading_AK4] > n_daus
    n_daughters_s_AK4 = vbf_data['Events']['L1PuppiAK4Jets_nDau'].array()[nonzero_AK4][subleading_AK4] > -1


    dau_mask_AK4 = n_daughters_l_AK4*n_daughters_s_AK4

    comb_gen_mask = dau_mask_gen
    comb_AK4_mask = dau_mask_AK4