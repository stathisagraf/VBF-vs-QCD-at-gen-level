import numpy as np
import matplotlib.pyplot as plt
import uproot 
import mplhep as mpl


def hist(ax,data,len_data, bins, density, color, label,y_scale_log = False, x_scale_log = False ):
    
   
    
    
    h = ax.hist(data,bins, histtype='step',density = density,color= color, label=label)
    bin_centers = 0.5 * (h[1][:-1] + h[1][1:])
    bin_width = h[1][1]-h[1][0]
    ax.legend()
    
    
    if y_scale_log ==True:
        ax.set_yscale('log')
    if x_scale_log == True:
        ax.set_xscale('log')
    
    
    for i in range(len_data):
        if density ==True:
            plt.errorbar(bin_centers, h[0][i], yerr=np.sqrt(h[0][i]/(bin_width*len(data[i]))), fmt='none', linewidth=1, color=color[i],capsize=2)
            
        if density==False:
            plt.errorbar(bin_centers, h[0][i], yerr=np.sqrt(h[0][i]), fmt='none', linewidth=1, color=color[i],capsize=2)
    return h  


