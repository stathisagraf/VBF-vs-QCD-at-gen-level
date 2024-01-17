import numpy as np
import matplotlib.pyplot as plt
import uproot 
import mplhep as mpl


def hist(ax,data,len_data, bins, density, color, label,y_scale_log = False, x_scale_log = False, legend_font = 11 ):
    
   
    
    
    h = ax.hist(data,bins, histtype='step',density = density,color= color, label=label)
    bin_centers = 0.5 * (h[1][:-1] + h[1][1:])
    bin_width = h[1][1]-h[1][0]
    ax.legend(fontsize = legend_font)
    
    
    if y_scale_log ==True:
        ax.set_yscale('log')
    if x_scale_log == True:
        ax.set_xscale('log')
    
    if len_data!=1:
        
        for i in range(len_data):
            if density ==True:
                plt.errorbar(bin_centers, h[0][i], yerr=np.sqrt(h[0][i]/(bin_width*len(data[i]))), fmt='none', linewidth=1, color=color[i],capsize=2)
                
            if density==False:
                plt.errorbar(bin_centers, h[0][i], yerr=np.sqrt(h[0][i]), fmt='none', linewidth=1, color=color[i],capsize=2)
    if len_data == 1:
        if density ==True:
            plt.errorbar(bin_centers, h[0], yerr=np.sqrt(h[0]/(bin_width*len(data))), fmt='none', linewidth=1, color=color,capsize=2)
            
        if density==False:
            plt.errorbar(bin_centers, h[0], yerr=np.sqrt(h[0]), fmt='none', linewidth=1, color=color,capsize=2)
        
    return h  


