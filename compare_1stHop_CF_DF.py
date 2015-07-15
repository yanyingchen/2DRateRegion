'''
Created on Apr 5, 2015

@author: yanying
'''

import coding_scheme as cs 

 

# for plotting 
import rate_region as rr
import matplotlib.pyplot as plt
#import numpy as np
import matplotlib.pylab as pylab

# for saving the figure 
import savefig as savefig 


class Data():
    '''
    keep track of the data computed
    '''
    def __init__(self, Ps, 
                g13, g14, g23, g24,
                font={'family' : 'serif',
                                   'color'  : 'darkred',
                                   'weight' : 'normal',
                                   'size'   : 16,}
                ):
        self.Ps = Ps
        self.g13 = g13
        self.g14 = g14 
        self.g23 = g23 
        self.g24 = g24
        self.font = font 
        
        self.df_Region = cs.DF_scheme(Ps, Ps, g13, g14, g23, g24) 
        self.cf_Region = cs.CF_scheme_integerCoeff(Ps, g13, g14, g23, g24)

def display(oneSimulation,pathString ,saveOnly=True):
    '''
    show/save figures
    '''
    self = oneSimulation

    pylab.rc('axes', linewidth=2) # make the axes boundary lines bold 
    #fig, ax = plt.subplots()
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.7, 0.7]) # left, bottom, width, height (range 0 to 1)
    rr.plot( self.df_Region, 'red', lw=6, axes=ax, label='DF')
    rr.plot( self.cf_Region, 'blue', lw=4, axes=ax, label='CF')  
    
    
    fig.suptitle('$g ={},{},{},{}$'.format(self.g13, self.g14, self.g23, self.g24), fontsize=14, fontweight='bold')
    ax.set_title(r'$P_s=P_1=P_2={}, \, N=1$'.format(self.Ps),fontdict=self.font)
    
    ax.set_xlabel('$R_1$', fontdict=self.font)
    ax.set_ylabel('$R_2$', fontdict=self.font)
    #ax.set_xlim(xmin=0, xmax=3.5) 
    #ax.set_ylim(ymin=0, ymax=3.5) 
    ax.legend(loc=0)
    
    nameStr = 'Ps_{}_g13g14g23g24_{}_{}_{}_{}'.format(self.Ps, self.g13, self.g14, self.g23, self.g24)
    savefig.save(path='{}/plots_DF_CF/{}'.format(pathString, nameStr), ext='pdf', close=saveOnly, verbose=True)
        



