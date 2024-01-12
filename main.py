import numpy as np
import matplotlib.pyplot as plt
import uproot 
import awkward as ak
import vector
import mplhep as hep
import time 
start = time.time()


from hist_plotting import hist #
from pull_vector_calc import pv #pv takes 5 arguments and outputs 5
from pull_vector_calc_dau_mass import pv_dau_mass
from relative_pull_vector_calc_dau_mass import rel_pv_dau_mass


plt.style.use(hep.style.CMS)

vbf_path = 'perfNano_VBFHInv_PU200.l1ctlayer1.root'
qcd_path = 'perfNano_QCD_PU200.l1ctlayer1.root'

# dau_pt_cut_table = [0,1,2,3,4,5,6,7,8,9,10]

# leading_jet_pt_cuts = [20,30,40,50]
# for pt_cut_leading in leading_jet_pt_cuts:
        
#     for dau_cuts_pt in dau_pt_cut_table:
            
#         vbf_pv_y_leading, vbf_pv_phi_leading, vbf_pv_y_subleading, vbf_pv_phi_subleading, x = pv(vbf_path,pt_cut_leading,0,0,0, dau_pt_cuts = dau_cuts_pt)
        
#         qcd_pv_y_leading, qcd_pv_phi_leading, qcd_pv_y_subleading, qcd_pv_phi_subleading, z = pv(qcd_path,pt_cut_leading,0,0,0, dau_pt_cuts = dau_cuts_pt)
        
        
#         vbf_pv_leading = vector.array({
#                                         "x": vbf_pv_y_leading,
#                                         "y": vbf_pv_phi_leading
#                                         })
        
#         vbf_pa_leading = vbf_pv_leading.phi
#         vbf_pv_magnitude_leading = vbf_pv_leading.rho
        
        
#         vbf_pv_subleading = vector.array({
#                                         "x": vbf_pv_y_subleading,
#                                         "y": vbf_pv_phi_subleading
#                                         })
        
#         vbf_pa_subleading = vbf_pv_subleading.phi
#         vbf_pv_magnitude_subleading = vbf_pv_subleading.rho
        
        
#         qcd_pv_leading = vector.array({
#                                         "x": qcd_pv_y_leading,
#                                         "y": qcd_pv_phi_leading
#                                         })
        
#         qcd_pa_leading = qcd_pv_leading.phi
#         qcd_pv_magnitude_leading = qcd_pv_leading.rho
        
        
#         qcd_pv_subleading = vector.array({
#                                         "x": qcd_pv_y_subleading,
#                                         "y": qcd_pv_phi_subleading
#                                         })
        
#         qcd_pa_subleading = qcd_pv_subleading.phi
#         qcd_pv_magnitude_subleading = qcd_pv_subleading.rho
        
        
#         fig,ax = plt.subplots()
        
        
#         hist(ax,[vbf_pa_leading,qcd_pa_leading],2,10,True,['black','red'],['VBF pull angle leading','QCD pull angle leading'] ,x_scale_log = False, y_scale_log = False )
#         plt.text(-1.8,0.1,r'Leading jet $p_T$ > {} GeV'.format(pt_cut_leading))
#         plt.text(-1.8,0.08,r'Daughter $p_T$ > {} GeV'.format(dau_cuts_pt))
#         plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Leading_Pull_angle_pt_cut_{}_dau_cuts_{}.pdf'.format(pt_cut_leading,dau_cuts_pt))
#         plt.show()
        
        
        
#         fig,ax = plt.subplots()
        
        
#         hist(ax,[vbf_pa_subleading,qcd_pa_subleading],2,10,True,['black','red'],['VBF pull angle subleading','QCD pull angle subleading'] ,x_scale_log = False, y_scale_log = False )
#         plt.text(-1.8,0.1,r'Leading jet $p_T$ > {} GeV'.format(pt_cut_leading))
#         plt.text(-1.8,0.08,r'Daughter $p_T$ > {} GeV'.format(dau_cuts_pt))
#         plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Subeading_Pull_angle_pt_cut_{}_dau_cuts_{}.pdf'.format(pt_cut_leading,dau_cuts_pt))
    
#         plt.show()
        

    
# dau_mass_cut_table = [0,0.1,0.11,0.15,0.5]
# dau_pt_cuts = [0,1,4,6,10,15]
# leading_jet_pt_cuts = [0,20,50,70]
# for pt_cut_leading in leading_jet_pt_cuts:
#     for dau_pt_cut in dau_pt_cuts:
#         for dau_cuts_mass in dau_mass_cut_table:
                
#             vbf_pv_y_leading, vbf_pv_phi_leading, vbf_pv_y_subleading, vbf_pv_phi_subleading, x = pv_dau_mass(vbf_path,pt_cut_leading,0,0,0,dau_pt_cuts = dau_pt_cut,  dau_mass_cuts = dau_cuts_mass)
            
#             qcd_pv_y_leading, qcd_pv_phi_leading, qcd_pv_y_subleading, qcd_pv_phi_subleading, z = pv_dau_mass(qcd_path,pt_cut_leading,0,0,0,dau_pt_cuts = dau_pt_cut,  dau_mass_cuts = dau_cuts_mass)
            
            
#             vbf_pv_leading = vector.array({
#                                             "x": vbf_pv_y_leading,
#                                             "y": vbf_pv_phi_leading
#                                             })
            
#             vbf_pa_leading = vbf_pv_leading.phi
#             vbf_pv_magnitude_leading = vbf_pv_leading.rho
            
            
#             vbf_pv_subleading = vector.array({
#                                             "x": vbf_pv_y_subleading,
#                                             "y": vbf_pv_phi_subleading
#                                             })
            
#             vbf_pa_subleading = vbf_pv_subleading.phi
#             vbf_pv_magnitude_subleading = vbf_pv_subleading.rho
            
            
#             qcd_pv_leading = vector.array({
#                                             "x": qcd_pv_y_leading,
#                                             "y": qcd_pv_phi_leading
#                                             })
            
#             qcd_pa_leading = qcd_pv_leading.phi
#             qcd_pv_magnitude_leading = qcd_pv_leading.rho
        
#             qcd_pv_subleading = vector.array({
#                                             "x": qcd_pv_y_subleading,
#                                             "y": qcd_pv_phi_subleading
#                                             })
            
#             qcd_pa_subleading = qcd_pv_subleading.phi
#             qcd_pv_magnitude_subleading = qcd_pv_subleading.rho
            
            
#             fig,ax = plt.subplots()
            
#             print(len(vbf_pa_leading),len(qcd_pa_leading))
#             hist(ax,[vbf_pa_leading,qcd_pa_leading],2,9,True,['black','red'],['VBF pull angle leading','QCD pull angle leading'] ,x_scale_log = False, y_scale_log = False )
#             plt.text(-1.8,0.1,r'Leading jet $p_T$ > {} GeV'.format(pt_cut_leading))
#             plt.text(-1.8,0.085,r'Daughter $mass$ > {} GeV'.format(dau_cuts_mass))
#             plt.text(-1.8,0.07,r'Daughter $p_T$ > {} GeV'.format(dau_pt_cut))
#             plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Leading_Pull_angle_pt_cut_{}_mass_cut_{}_dau_cuts_{}.pdf'.format(dau_pt_cut,pt_cut_leading,dau_cuts_mass))
#             plt.show()
            
            
            
#             fig,ax = plt.subplots()
            
            
#             hist(ax,[vbf_pa_subleading,qcd_pa_subleading],2,9,True,['black','red'],['VBF pull angle subleading','QCD pull angle subleading'] ,x_scale_log = False, y_scale_log = False )
#             plt.text(-1.8,0.1,r'Leading jet $p_T$ > {} GeV'.format(pt_cut_leading))
#             plt.text(-1.8,0.085,r'Daughter $mass$ > {} GeV'.format(dau_cuts_mass))
#             plt.text(-1.8,0.07,r'Daughter $p_T$ > {} GeV'.format(dau_pt_cut))
#             plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Subeading_Pull_angle_pt_cut_{}_mass_cut_{}_dau_cuts_{}.pdf'.format(dau_pt_cut,pt_cut_leading,dau_cuts_mass))
        
#             plt.show()
    
    


dau_mass_cut_table = [0,0.1,0.11,0.15,0.5]
dau_pt_cuts = [0,1,4,6,10,15]
leading_jet_pt_cuts = [0,20,50,70]
for pt_cut_leading in leading_jet_pt_cuts:
    for dau_pt_cut in dau_pt_cuts:
        for dau_cuts_mass in dau_mass_cut_table:
                
            vbf_rel_pa_leading, vbf_rel_pa_subleading, x = rel_pv_dau_mass(vbf_path,pt_cut_leading,0,0,0,dau_pt_cuts = dau_pt_cut,  dau_mass_cuts = dau_cuts_mass)
            
            qcd_rel_pa_leading, qcd_rel_pa_subleading, z = rel_pv_dau_mass(qcd_path,pt_cut_leading,0,0,0,dau_pt_cuts = dau_pt_cut,  dau_mass_cuts = dau_cuts_mass)
            
            
            
            
            fig,ax = plt.subplots()
            
            print(len(vbf_rel_pa_leading),len(qcd_rel_pa_leading))
            hist(ax,[vbf_rel_pa_leading,qcd_rel_pa_leading],2,9,True,['black','red'],['VBF relative pull angle leading','QCD relative pull angle leading'] ,x_scale_log = False, y_scale_log = False )
            plt.text(-1.8,0.1,r'Leading jet $p_T$ > {} GeV'.format(pt_cut_leading))
            plt.text(-1.8,0.085,r'Daughter $mass$ > {} GeV'.format(dau_cuts_mass))
            plt.text(-1.8,0.07,r'Daughter $p_T$ > {} GeV'.format(dau_pt_cut))
            plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Leading_Rel_pull_angle_pt_cut_{}_mass_cut_{}_dau_cuts_{}.pdf'.format(dau_pt_cut,pt_cut_leading,dau_cuts_mass))
            plt.show()
            
            
            
            fig,ax = plt.subplots()
            
            
            hist(ax,[vbf_rel_pa_subleading,qcd_rel_pa_subleading],2,9,True,['black','red'],['VBF relative pull angle subleading','QCD relative pull angle subleading'] ,x_scale_log = False, y_scale_log = False )
            plt.text(-1.8,0.1,r'Leading jet $p_T$ > {} GeV'.format(pt_cut_leading))
            plt.text(-1.8,0.085,r'Daughter $mass$ > {} GeV'.format(dau_cuts_mass))
            plt.text(-1.8,0.07,r'Daughter $p_T$ > {} GeV'.format(dau_pt_cut))
            plt.savefig('/home/stathis/Desktop/Research project/Data/VBF vs QCD at gen/Plots/Subeading_Rel_pull_angle_pt_cut_{}_mass_cut_{}_dau_cuts_{}.pdf'.format(dau_pt_cut,pt_cut_leading,dau_cuts_mass))
        
            plt.show()










end = time.time()

print(end-start, ' seconds')