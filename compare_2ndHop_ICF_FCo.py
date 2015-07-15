'''
Created on Apr 6, 2015

@author: yanying
'''


import coding_scheme as cs 
import alternative_compute_ICFsch3 as alternative
 

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
    def __init__(self, P3, P4, 
                numPoints,
                font={'family' : 'serif',
                                   'color'  : 'darkred',
                                   'weight' : 'normal',
                                   'size'   : 16,}
                ):
        self.P3 = P3
        self.P4 = P4
        self.font = font 
        
        # self.icf_Region_bigR1 = alternative.ICF_sch3_bigR1(P3, P4, numPoints)
        self.icf_Region_bigR2 = alternative.ICF_sch3_bigR2(P3, P4, numPoints)
        self.icf_Region_bigR1 = cs._ICF_sch3_bigR2_fromBigR1Region(self.icf_Region_bigR2)
        
        self.fco_Region = cs.FCo_scheme(P3, P4)

def display(oneSimulation,pathString ,saveOnly=True):
    '''
    show/save figures
    '''
    self = oneSimulation

    pylab.rc('axes', linewidth=2) # make the axes boundary lines bold 
    #fig, ax = plt.subplots()
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.7, 0.7]) # left, bottom, width, height (range 0 to 1)
    rr.plot( self.fco_Region, 'red', lw=6, axes=ax, label='FCo')
    rr.plot( self.icf_Region_bigR1, 'blue', lw=4, axes=ax, label='ICF')  
    rr.plot( self.icf_Region_bigR2, 'blue', lw=4, axes=ax)  
        
    ax.set_title(r'$P_3={}, \, P_4={}, \, N=1$'.format(self.P3, self.P4),fontdict=self.font)
    
    ax.set_xlabel('$R_1$', fontdict=self.font)
    ax.set_ylabel('$R_2$', fontdict=self.font)
    #ax.set_xlim(xmin=0, xmax=3.5) 
    #ax.set_ylim(ymin=0, ymax=3.5) 
    ax.legend(loc=0)
    
    nameStr = 'P3P4_{}_{}'.format(self.P3, self.P4)
    savefig.save(path='{}/plots_ICF_FCo/{}'.format(pathString, nameStr), ext='pdf', close=saveOnly, verbose=True)
        



