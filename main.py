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




plt.style.use(hep.style.CMS)

vbf_path = 'perfNano_VBFHInv_PU200.l1ctlayer1.root'
qcd_path = 'perfNano_QCD_PU200.l1ctlayer1.root'


vbf_pv_y_leading, vbf_pv_phi_leading, vbf_pv_y_subleading, vbf_pv_phi_subleading, x = pv(vbf_path,20,0,0,0)

qcd_pv_y_leading, qcd_pv_phi_leading, qcd_pv_y_subleading, qcd_pv_phi_subleading, z = pv(qcd_path,20,0,0,0)


vbf_pv_leading = vector.array({
                                "x": vbf_pv_y_leading,
                                "y": vbf_pv_phi_leading
                                })

vbf_pa_leading = vbf_pv_leading.phi
vbf_pv_magnitude_leading = vbf_pv_leading.rho


vbf_pv_subleading = vector.array({
                                "x": vbf_pv_y_subleading,
                                "y": vbf_pv_phi_subleading
                                })

vbf_pa_subleading = vbf_pv_subleading.phi
vbf_pv_magnitude_subleading = vbf_pv_subleading.rho


qcd_pv_leading = vector.array({
                                "x": qcd_pv_y_leading,
                                "y": qcd_pv_phi_leading
                                })

qcd_pa_leading = qcd_pv_leading.phi
qcd_pv_magnitude_leading = qcd_pv_leading.rho


qcd_pv_subleading = vector.array({
                                "x": qcd_pv_y_subleading,
                                "y": qcd_pv_phi_subleading
                                })

qcd_pa_subleading = qcd_pv_subleading.phi
qcd_pv_magnitude_subleading = qcd_pv_subleading.rho


fig,ax = plt.subplots()


hist(ax,[vbf_pa_leading,qcd_pa_leading],2,100,True,['black','red'],['VBF pull angle leading','QCD pull angle leading'] ,x_scale_log = False, y_scale_log = False )

plt.show()






















end = time.time()

print(end-start, ' seconds')