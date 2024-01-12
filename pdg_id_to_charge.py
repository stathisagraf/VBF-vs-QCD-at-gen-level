import numpy as np
import uproot 
import awkward as ak
import matplotlib.pyplot as plt
from collections import Counter

def charge(charg):
    return 0
    
    
vbf_path = 'perfNano_VBFHInv_PU200.l1ctlayer1.root'

data = uproot.open(vbf_path)
cha = []
for dau in range(5):

    cha.append(ak.flatten(data['Events']['GenJets_dau{}_pdgId'.format(dau)].array()))


count = Counter((ak.flatten(cha)))

print(count)
