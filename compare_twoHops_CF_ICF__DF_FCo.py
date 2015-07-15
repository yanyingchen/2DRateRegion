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
    def __init__(self, Ps, 
                 P3, P4, 
                numPoints,
                g13, g14, g23, g24,
                font={'family' : 'serif',
                                   'color'  : 'darkred',
                                   'weight' : 'normal',
                                   'size'   : 16,}
                ):
        self.Ps = Ps
        self.P3 = P3
        self.P4 = P4
        self.g13 = g13
        self.g14 = g14 
        self.g23 = g23 
        self.g24 = g24
        self.font = font 
        
        self.df_Region = cs.DF_scheme(Ps, Ps, g13, g14, g23, g24) 
        self.cf_Region = cs.CF_scheme_integerCoeff(Ps, g13, g14, g23, g24)
        
        self.icf_Region_bigR2 = alternative.ICF_sch3_bigR2(P3, P4, numPoints)
        self.icf_Region_bigR1 = cs._ICF_sch3_bigR2_fromBigR1Region(self.icf_Region_bigR2)
        
        self.fco_Region = cs.FCo_scheme(P3, P4)
        
        self.cf_icf_region = cs.TwoHops_CF_ICFsch3_fromRegions(self.cf_Region, 
                                                               self.icf_Region_bigR1, self.icf_Region_bigR2)
        self.df_fco_region = cs.TwoHops_DFIntegerCoeff_FCo_fromRegions(self.df_Region, self.fco_Region)
        
        
def display(oneSimulation,pathString, separate=False ,saveOnly=True):
    '''
    show/save figures
    '''
    self = oneSimulation

    pylab.rc('axes', linewidth=2) # make the axes boundary lines bold 
    #fig, ax = plt.subplots()
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.7, 0.7]) # left, bottom, width, height (range 0 to 1)
    if separate:
        rr.plot( self.df_Region, 'red', lw=6, axes=ax, label='DF')
        rr.plot( self.cf_Region, 'blue', lw=4, axes=ax, label='CF')  
        rr.plot( self.fco_Region, 'red', lw=6, axes=ax, label='DF')
        rr.plot( self.icf_Region_bigR1, 'blue', lw=4, axes=ax, label='ICF')  
        rr.plot( self.icf_Region_bigR2, 'blue', lw=4, axes=ax)
    
    rr.plot(self.df_fco_region, 'red', lw=6, label='DF+FCo')
    rr.plot(self.cf_icf_region, 'blue', lw=4, label='CF+ICF')
    fig.suptitle('$g ={},{},{},{}$'.format(self.g13, self.g14, self.g23, self.g24), fontsize=14, fontweight='bold')
    ax.set_title(r'$P_s=P_1=P_2={}, \, P_3={}, \, P_4={}, \, N=1$'.format(self.Ps,
                 self.P3, self.P4),fontdict=self.font)
    
    ax.set_xlabel('$R_1$', fontdict=self.font)
    ax.set_ylabel('$R_2$', fontdict=self.font)
    #ax.set_xlim(xmin=0, xmax=3.5) 
    #ax.set_ylim(ymin=0, ymax=3.5) 
    ax.legend(loc=0)
    
    nameStr = 'PsP3P4_{}_{}_{}_g13g14g23g24_{}_{}_{}_{}'.format(self.Ps, self.P3, self.P4, 
                                                                self.g13, self.g14, self.g23, self.g24)
    savefig.save(path='{}/plots_CF_ICF__DF_FCo/{}'.format(pathString, nameStr), ext='pdf', close=saveOnly, verbose=True)        
        