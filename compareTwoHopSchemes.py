'''
Created on Mar 31, 2015

compare three two-hop schemes
CF + ICFsch3
DF + FCo
p2p + MAC
@author: yanying
'''


import coding_scheme as cs 
import rate_region as rr 

# for plotting 
import matplotlib.pyplot as plt
#import numpy as np
import matplotlib.pylab as pylab
# from textwrap import wrap 


# for saving the figure 
import savefig as savefig 

class TwoHopSchemes(object):
    '''
    compute the rate region for the network by three TwoHop schemes.
    '''
    def __init__(self, Ps, P3, P4,
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
        
        self.region_CF_ICFsch3 = cs.TwoHops_CF_ICFsch3(Ps, P3, P4, g13, g14, g23, g24)
        self.region_DFIntegerCoeff_FCo = cs.TwoHops_DFIntegerCoeff_FCo(Ps, P3, P4, g13, g14, g23, g24)
        self.region_noInterference = cs.TwoHops_noInterference(Ps, Ps, P3, P4)
        
        
def display_TwoHopSchemes(two_hop_schemes, pathString):
    '''
    Produce a graph that compare ... 
    '''
    self = two_hop_schemes

    pylab.rc('axes', linewidth=2) # make the axes boundary lines bold 
    #fig, ax = plt.subplots()
    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.7, 0.7]) # left, bottom, width, height (range 0 to 1)
    rr.plot( self.region_CF_ICFsch3, 'red', lw=8, axes=ax, label='CF + ICF')
    rr.plot( self.region_DFIntegerCoeff_FCo, 'c', lw=4, axes=ax, label='DF + FCo')
    rr.plot( self.region_noInterference, 'blue', lw=2, axes=ax, label='Free Interf.')  
    
    
    fig.suptitle('$g ={},{},{},{}$'.format(self.g13, self.g14, self.g23, self.g24), fontsize=14, fontweight='bold')
    ax.set_title(r'$P_s=P_1=P_2={}, \, P_3={}, \, P_4={}, \, N=1$'.format(self.Ps, self.P3, self.P4),fontdict=self.font)
    
    ax.set_xlabel('$R_1$', fontdict=self.font)
    ax.set_ylabel('$R_2$', fontdict=self.font)
    #ax.set_xlim(xmin=0, xmax=3.5) 
    #ax.set_ylim(ymin=0, ymax=3.5) 
    ax.legend(loc=0)
    
    # savefig.save(path='{}compare_TwoHop_Schemes/PsP3P4_{}_{}_{}_g_{}_{}_{}_{}'.format(pathString, self.Ps, self.P3, self.P4, self.g13, self.g14, self.g23, self.g24), ext='pdf', close=False, verbose=True)    
    
    nameStr = 'PsP3P4_{}_{}_{}_g13g14g23g24_{}_{}_{}_{}'.format(self.Ps, self.P3, self.P4, self.g13, self.g14, self.g23, self.g24)
    savefig.save(path='{}compare_TwoHop_Schemes/{}'.format(pathString, nameStr), ext='pdf', close=False, verbose=True)
    
        
        