'''
Created on Apr 5, 2015

Returns the region (a numpy array) for ICF Scheme 3.
@author: yanying
'''

import alternative_rate_region as alternative_rr 
import rate_region as rr 

import shapely.geometry as geo

import numpy as np 



def _computeC(x):
    '''
    compute the value of a C( ) function, which denotes the capacity of a 
    point-to-point Gaussian channel
    '''
    return 0.5*np.log2(1+x)

def _ICF_sch3_computation(P3, P4, b1, b2):
    '''
    A coherent scheme with cardinality bounding.
    Compute the right hand side values
    ''' 
    r_min = np.min([_computeC( (1-b1)*P3 ),
                    _computeC( (1-b2)*P4 ),
                    0.5*_computeC( (1-b1)*P3+(1-b2)*P4 )
                    ]) 
    r_sum = _computeC( P3+P4+2*np.sqrt(b1*b2*P3*P4) )
    return [r_min, r_sum]


def _ICF_sch3_bigR1_onePowerAllocation(P3, P4, b1, b2):
    '''
    A coherent scheme with cardinality bounding.
    Note: Rmin is R2
    '''
    r_min, r_sum = _ICF_sch3_computation(P3, P4, b1, b2)
    inequalities = [(-1, 1, 0), # -R1 + R2 <= 0  
                    (0, 1, r_min), #  0 + R2  <= r_min
                    (1, 1, r_sum)] # R1 + R2 <= r_sum 
#     return rr.InequalityRegion(inequalities)
    return alternative_rr.ICFsch3RateRegion(inequalities)
    

def ICF_sch3_bigR1(P3, P4, numPoints=21):
    '''
    Returns 
    
    A coherent scheme with cardinality bounding.
    Note: Rmin is R2
    '''
    rateRegionList = [_ICF_sch3_bigR1_onePowerAllocation(P3, P4, x, y) for x in np.linspace(0, 1, numPoints) for y in np.linspace(0, 1, numPoints) ]
    points =  alternative_rr.unionAndConvexhull(rateRegionList)
    return rr.RateRegion( geo.MultiPoint(points).convex_hull)

def _ICF_sch3_bigR2_onePowerAllocation(P3, P4, b1, b2):
    '''
    A coherent scheme with cardinality bounding.
    Note: Rmin is R1
    '''
    r_min, r_sum = _ICF_sch3_computation(P3, P4, b1, b2)
    inequalities = [(1, -1, 0), # R1 - R2 <= 0 ...........   
                    (1, 0, r_min), #  R1 + 0  <= r_min
                    (1, 1, r_sum)] # R1 + R2 <= r_sum 
#     return rr.InequalityRegion(inequalities) 
    return alternative_rr.ICFsch3RateRegion(inequalities)

def ICF_sch3_bigR2(P3, P4, numPoints=21):
    '''
    A coherent scheme with cardinality bounding.
    Note: Rmin is R2
    '''            
    rateRegionList = [_ICF_sch3_bigR2_onePowerAllocation(P3, P4, x, y) for x in np.linspace(0, 1, numPoints) for y in np.linspace(0, 1, numPoints) ]
    points =  alternative_rr.unionAndConvexhull(rateRegionList)
    return rr.RateRegion( geo.MultiPoint(points).convex_hull)