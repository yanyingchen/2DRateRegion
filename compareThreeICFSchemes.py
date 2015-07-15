'''
Created on Mar 30, 2015

@author: yanying
'''

import coding_scheme as cs 
import rate_region as rr 

# for plotting 
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pylab as pylab

# for saving the figure 
import savefig as savefig 

class Data(object):
    '''
    compare three ICF schemes
    '''
    def __init__(self, P3, P4,font={'family' : 'serif',
                                    'color'  : 'darkred',
                                    'weight' : 'normal',
                                    'size'   : 16,}
                 ):
        self.P3 = P3
        self.P4 = P4
        self.font = font 
        
        # get the rate regions for bigR1
        self.icf_sch1_bigR1 = cs.ICF_sch1_bigR1(P3, P4)
        self.icf_sch2_bigR1 = cs.ICF_sch2_bigR1(P3, P4)
        self.icf_sch3_bigR1 = cs.ICF_sch3_bigR1(P3, P4)
        
        # get the rate regions for bigR2
        self.icf_sch1_bigR2 = self.icf_sch1_bigR1.symmetricRegion_yEqualToX()
        self.icf_sch2_bigR2 = self.icf_sch2_bigR1.symmetricRegion_yEqualToX()
        self.icf_sch3_bigR2 = self.icf_sch3_bigR1.symmetricRegion_yEqualToX()


def display(oneSimulation,pathString ,saveOnly=True):
    '''
    Produce a graph that compare three ICF schemes
    '''
    self = oneSimulation
    
    # get the union rate region for plotting 
    icf_sch1 = rr.union( [self.icf_sch1_bigR1, self.icf_sch1_bigR2] )
    icf_sch2 = rr.union( [self.icf_sch2_bigR1, self.icf_sch2_bigR2] )
    icf_sch3 = rr.union( [self.icf_sch3_bigR1, self.icf_sch3_bigR2] )

    pylab.rc('axes', linewidth=2) # make the axes boundary lines bold 
    fig, ax = plt.subplots()
    rr.plot( icf_sch1, 'g', axes=ax, label='Scheme 1')
    rr.plot( icf_sch2, 'b', axes=ax, label='Scheme 2')
    rr.plot( icf_sch3, 'r', axes=ax, label='Scheme 3')  
    
    # plot the line: R2 = R1
    tmp = np.asarray(self.icf_sch1_bigR1._geometry.boundary)
    tmp2 = [[0, tmp[1, 0]], [0, tmp[1, 1] ] ]
    ax.plot( [0, tmp[1, 0]], [0, tmp[1, 1] ] , 'k--', lw=2)

    ax.set_title(r'$ P_3={}, \, P_4={}, \, N=1$'.format(self.P3, self.P4) , fontdict=self.font)
    ax.set_xlabel('$R_1$', fontdict=self.font)
    ax.set_ylabel('$R_2$', fontdict=self.font)
    ax.set_xlim(xmin=0, xmax=3.5) 
    ax.set_ylim(ymin=0, ymax=3.5) 
    ax.legend(loc='upper right')
    
    savefig.save(path='{}/compare_three_ICF_Schemes/P3P4_{}_{}'.format(pathString, self.P3, self.P4 ), ext='pdf', close=saveOnly, verbose=True)
    if (0):
        # add annotations for 3 regions
        plt.text(0.7, 2.3, 'capacity region by coherent coding with cardinality-bounding', color='red')
        plt.text(1, 1.8, 'non-coherent coding with cardinality-bounding', color='blue')
        plt.text(1.3, 1.4, 'capacity region by coherent coding with cardinality-bounding', color='green')
        
        bbox_props = dict(boxstyle="round,pad=0.1", fc="white", ec="g", lw=1)
        plt.text(1, 0.5, r'$1$', color='black', bbox=bbox_props)
        bbox_props = dict(boxstyle="round,pad=0.1", fc="white", ec="b", lw=1)
        plt.text(1.7, 0.4, r'$2$', color='black', bbox=bbox_props)
        bbox_props = dict(boxstyle="round,pad=0.1", fc="white", ec="r", lw=1)
        plt.text(2.2, 0.3, r'$3$', color='black', bbox=bbox_props)
    
    # plot the line: R2 = R1
    # plt.plot( [0, xyArray1[1, 0]], [0, xyArray1[1, 1] ] , 'k--', lw=2) 

    # return pFig
