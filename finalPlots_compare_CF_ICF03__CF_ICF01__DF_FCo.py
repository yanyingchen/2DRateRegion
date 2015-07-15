'''
Created on Apr 6, 2015

@author: yanying
'''


import coding_scheme as cs 
import alternative_compute_ICFsch3 as alternative
 

# for plotting 
import rate_region as rr
import matplotlib.pyplot as plt
import numpy as np
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
        
        # first hop:
        self.df_region = cs.DF_scheme(Ps, Ps, g13, g14, g23, g24) 
        self.cf_region = cs.CF_scheme_integerCoeff(Ps, g13, g14, g23, g24)
        # second hop:
        self.icf_region_bigR2 = alternative.ICF_sch3_bigR2(P3, P4, numPoints) # ICF Scheme 3
        self.icf_region_bigR1 = cs._ICF_sch3_bigR2_fromBigR1Region(self.icf_region_bigR2)        
        self.fco_region = cs.FCo_scheme(P3, P4) # FCo
        self.icf01_region_bigR1 = cs.ICF_sch1_bigR1(P3, P4) # ICF Scheme 1
        self.icf01_region_bigR2 = cs.ICF_sch1_bigR2(P3, P4)
        # two hops:
        self.cf_icf_region = cs.TwoHops_CF_ICFsch3_fromRegions(self.cf_region, 
                                                               self.icf_region_bigR1, self.icf_region_bigR2)
        self.df_fco_region = cs.TwoHops_DFIntegerCoeff_FCo_fromRegions(self.df_region, self.fco_region)
        self.cf_icf01_region = cs.TwoHops_CF_ICFsch1_fromRegions(self.icf01_region_bigR1, self.icf01_region_bigR2, self.cf_region) 
        
        
def display70(oneSimulation,pathString, separate=False ,saveOnly=True):
    '''
    show/save figures
    '''
    self = oneSimulation
    
    tmp = np.asarray(oneSimulation.cf_icf01_region._geometry.boundary)
    tmp01= tmp[(0,1,3), :]
    
    tmp = np.asarray(oneSimulation.cf_icf_region._geometry.boundary)
    tmp03= tmp[(1,2,3,4), :]
    
    tmp = np.asarray(oneSimulation.df_fco_region._geometry.boundary)
    tmpDF = tmp[(0,1,2,3), :]


    

    pylab.rc('axes', linewidth=2) # make the axes boundary lines bold 
#     fig, ax = plt.subplots()
    fig = plt.figure()
    fig.set( size_inches=(8.8, 6) )
    
    ax = fig.add_axes([0.1, 0.1, 0.7, 0.7]) # left, bottom, width, height (range 0 to 1)
#     if separate:
#         rr.plot( self.df_region, 'red', lw=6, axes=ax, label='DF')
#         rr.plot( self.cf_region, 'blue', lw=4, axes=ax, label='CF')  
#         rr.plot( self.fco_region, 'red', lw=6, axes=ax, label='DF')
#         rr.plot( self.icf_region_bigR1, 'blue', lw=4, axes=ax, label='ICF')  
#         rr.plot( self.icf_region_bigR2, 'blue', lw=4, axes=ax)
    
    
    
#     rr.plot(self.df_fco_region, 'blue', lw=6, label='DF+FCo')
#     rr.plot(self.cf_icf_region, 'red', lw=4, label='CF+ICF Scheme 3')
#     rr.plot(self.cf_icf01_region, 'green', lw=2, label='CF+ICF Scheme 1' )
    
    
    ax.plot(tmpDF[:,0], tmpDF[:,1], 'blue', lw=6, label='DF+FCo')
    ax.plot(tmp03[:,0], tmp03[:,1], 'red', lw=4, label='CF+ICF Scheme 3')
    ax.plot(tmp01[:,0], tmp01[:,1], 'green', lw=2, label='CF+ICF Scheme 1' )
    
    
    
    
#     fig.suptitle('$g ={},{},{},{}$'.format(self.g13, self.g14, self.g23, self.g24), fontsize=14, fontweight='bold')
    ax.set_title(r'$P_s=P_1=P_2={}, \, P_3={}, \, P_4={}, \, N=1$'.format(self.Ps,
                 self.P3, self.P4),fontdict=self.font)
    
    ax.set_xlabel('$R_1$', fontdict=self.font)
    ax.set_ylabel('$R_2$', fontdict=self.font)
    #ax.set_xlim(xmin=0, xmax=3.5) 
    #ax.set_ylim(ymin=0, ymax=3.5) 
    ax.legend(loc=0)
    
    nameStr = 'PsP3P4_{}_{}_{}_g13g14g23g24_{}_{}_{}_{}'.format(self.Ps, self.P3, self.P4, 
                                                                self.g13, self.g14, self.g23, self.g24)
    savefig.save(path='{}/plots_CF_ICF03__CF_ICF01__DF_FCo/{}'.format(pathString, nameStr), ext='pdf', close=saveOnly, verbose=True)        

def display2000(oneSimulation,pathString, separate=False ,saveOnly=True):
    '''
    show/save figures
    '''
    self = oneSimulation
    
    tmp = np.asarray(oneSimulation.df_fco_region._geometry.boundary)
    tmpDF = tmp[(0,1,2,3), :]
    
    tmp = np.asarray(oneSimulation.cf_icf01_region._geometry.boundary)
    tmp01= tmp[(0,1,2), :]
    
    tmp = np.asarray(oneSimulation.cf_icf_region._geometry.boundary)
    tmp03= tmp[(1,2,3,4), :]


    

    pylab.rc('axes', linewidth=2) # make the axes boundary lines bold 
#     fig, ax = plt.subplots()
    fig = plt.figure()
    fig.set( size_inches=(8.8, 6) )
    
    ax = fig.add_axes([0.1, 0.1, 0.7, 0.7]) # left, bottom, width, height (range 0 to 1)
#     if separate:
#         rr.plot( self.df_region, 'red', lw=6, axes=ax, label='DF')
#         rr.plot( self.cf_region, 'blue', lw=4, axes=ax, label='CF')  
#         rr.plot( self.fco_region, 'red', lw=6, axes=ax, label='DF')
#         rr.plot( self.icf_region_bigR1, 'blue', lw=4, axes=ax, label='ICF')  
#         rr.plot( self.icf_region_bigR2, 'blue', lw=4, axes=ax)
    
    
    
#     rr.plot(self.df_fco_region, 'blue', lw=6, label='DF+FCo')
#     rr.plot(self.cf_icf_region, 'red', lw=4, label='CF+ICF Scheme 3')
#     rr.plot(self.cf_icf01_region, 'green', lw=2, label='CF+ICF Scheme 1' )
    
    
    ax.plot(tmpDF[:,0], tmpDF[:,1], 'blue', lw=6, label='DF+FCo')
    ax.plot(tmp03[:,0], tmp03[:,1], 'red', lw=4, label='CF+ICF Scheme 3')
    ax.plot(tmp01[:,0], tmp01[:,1], 'green', lw=2, label='CF+ICF Scheme 1' )
    
    
    
    
#     fig.suptitle('$g ={},{},{},{}$'.format(self.g13, self.g14, self.g23, self.g24), fontsize=14, fontweight='bold')
    ax.set_title(r'$P_s=P_1=P_2={}, \, P_3={}, \, P_4={}, \, N=1$'.format(self.Ps,
                 self.P3, self.P4),fontdict=self.font)
    
    ax.set_xlabel('$R_1$', fontdict=self.font)
    ax.set_ylabel('$R_2$', fontdict=self.font)
    #ax.set_xlim(xmin=0, xmax=3.5) 
    #ax.set_ylim(ymin=0, ymax=3.5) 
    ax.legend(loc=3)
    
    nameStr = 'PsP3P4_{}_{}_{}_g13g14g23g24_{}_{}_{}_{}'.format(self.Ps, self.P3, self.P4, 
                                                                self.g13, self.g14, self.g23, self.g24)
    savefig.save(path='{}/plots_CF_ICF03__CF_ICF01__DF_FCo/{}'.format(pathString, nameStr), ext='pdf', close=saveOnly, verbose=True)        
