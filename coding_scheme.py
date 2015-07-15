'''
Created on Mar 29, 2015

Functions that return the rate region for a channel by first computing a list of 
inequalities in their canonical format: x*c0+y*c1<=c2

Convention: 
Both function and rate region are named after the coding scheme. 
the function name starts with upper case letter
the rate region variable starts with lower case letter  
@author: yanying
'''

import rate_region as rr 
import alternative_compute_ICFsch3 as alternative_compute 

# import shapely.geometry as geo 

import numpy as np 

def _computeC(x):
    '''
    compute the value of a C( ) function, which denotes the capacity of a 
    point-to-point Gaussian channel
    '''
    return 0.5*np.log2(1+x)

    
def ICF_sch1_bigR1(P3, P4):
    '''
    A non-coherent scheme without cardinality bounding. 
    Note: Rmax is R1
    '''
    r1 = np.min([_computeC( P3), 
                 _computeC( P4 ), 
                 0.5*_computeC( P3+P4 ) ])
    inequalities = [(-1, 1, 0), # -R1 + R2 <= 0  
                    (1, 0, r1)] #  R1 + 0  <= r1
    return rr.InequalityRegion(inequalities)

        
def ICF_sch1_bigR2(P3, P4):
    '''
    A non-coherent scheme without cardinality bounding.  
    Note: Rmax is R2 
    '''
    r2 = np.min([_computeC( P3 ), 
                 _computeC( P4 ), 
                 0.5*_computeC( P3+P4 ) ])
    inequalities = [(1, -1, 0), # R1 - R2 <= 0  
                    (0, 1, r2)] # 0  + R2 <= r2
    return rr.InequalityRegion(inequalities)   
        

def ICF_sch2_bigR1(P3, P4):
    '''
    A non-coherent scheme with cardinality bounding.
    Note: Rmin is R2
    '''
    r2 = np.min([_computeC( P3), _computeC( P4 )  ])
    rSum = _computeC( P3+P4 )
    inequalities = [(-1, 1, 0), # -R1 + R2 <= 0  
                    (0, 1, r2), #  0 + R2  <= r2
                    (1, 1, rSum)] # R1 + R2 <= rSum 
    return rr.InequalityRegion(inequalities)    


def ICF_sch2_bigR2(P3, P4):
    '''
    A non-coherent scheme with cardinality bounding.
    Note: Rmin is R1
    '''
    r1 = np.min([_computeC( P3), _computeC( P4 )  ])
    rSum = _computeC( P3+P4 )
    inequalities = [(1, -1, 0), # R1 - R2 <= 0  
                    (1, 0, r1), # R1 + 0  <= r1
                    (1, 1, rSum)] # R1 + R2 <= rSum 
    return rr.InequalityRegion(inequalities) 


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
    return rr.InequalityRegion(inequalities)


def ICF_sch3_bigR1(P3, P4, numPoints=21):
    '''
    A coherent scheme with cardinality bounding.
    Note: Rmin is R2
    '''
    rateRegionList = [_ICF_sch3_bigR1_onePowerAllocation(P3, P4, x, y) for x in np.linspace(0, 1, numPoints) for y in np.linspace(0, 1, numPoints) ]
    return rr.unionAndConvexhull(rateRegionList)


def _ICF_sch3_bigR2_onePowerAllocation(P3, P4, b1, b2):
    '''
    A coherent scheme with cardinality bounding.
    Note: Rmin is R1
    '''
    r_min, r_sum = _ICF_sch3_computation(P3, P4, b1, b2)
    inequalities = [(1, -1, 0), # R1 - R2 <= 0 ...........   
                    (1, 0, r_min), #  R1 + 0  <= r_min
                    (1, 1, r_sum)] # R1 + R2 <= r_sum 
    return rr.InequalityRegion(inequalities)
        
def ICF_sch3_bigR2(P3, P4, numPoints=101):
    '''
    A coherent scheme with cardinality bounding.
    Note: Rmin is R2
    '''            
    rateRegionList = [_ICF_sch3_bigR2_onePowerAllocation(P3, P4, x, y) for x in np.linspace(0, 1, numPoints) for y in np.linspace(0, 1, numPoints) ]      
    return rr.unionAndConvexhull(rateRegionList)
    
    
def _ICF_sch3_bigR2_fromBigR1Region(icf_sch3_bigR1):
    '''
    By symmetricity, 
    '''    
    return icf_sch3_bigR1.symmetricRegion_yEqualToX()
    
    
def CF_scheme_integerCoeff(Ps, g13, g14, g23, g24):
    '''
    Bobark's Compute-and-Forward, when the channel coefficients are integers. 
    I.e., equations that relays decode are matched to the channel coefficients.
    Notation convention: 
    g13: from (source) node 1 to (relay) node 3.
    Note: when all channel coefficients are non-zero! 
    '''    
    # constraints from the first relay node (node 3)
    # R1, R2 <= 0.5*np.log2( inverse_of_(g13^2 + g23^2 ) +Ps)
    # constraints from the second relay node (node 4)
    # R1, R2 <= 0.5*np.log2( inverse_of_(g14^2 + g24^2 ) +Ps)
    r_node3 = max(0.5*np.log2( 1/(g13**2 + g23**2 ) +Ps), 0)
    r_node4 = max(0.5*np.log2( 1/(g14**2 + g24**2 ) +Ps), 0)
    if r_node3==0:
        print('r_node3 is:', r_node3)
    if r_node4==0:
        print('r_node3 is:', r_node4)
    r_max = min(r_node3, r_node4)
    inequalities = [(1, 0, r_max ), # R1 + 0 <= r_max
                    (0, 1, r_max)]  # 0 + R2 <= r_max
    return rr.InequalityRegion(inequalities)


def MAC_capacity_independentMessages(P1, P2, h1=1, h2=1):
    '''
    capacity region of sending two independent messages over a two-user MAC.
    '''
    r1 = _computeC( (h1**2)*P1 )
    r2 = _computeC( (h2**2)*P2 )
    rSum = _computeC( (h1**2)*P1+(h2**2)*P2 )
    inequalities = [(1, 0, r1),  # R1 + 0  <= r1
                    (0, 1, r2),  # 0  + R2 <= r2
                    (1, 1, rSum)] # R1 + R2 <= rSum
    return  rr.InequalityRegion(inequalities)
   
    
def DF_scheme(P1, P2, g13, g14, g23, g24):
    '''
    Let two relay nodes (node 3 and node 4) both decode both messages, which 
    are independent.
    DF region is the intersection of two 2-user MAC capacity region
    For details: check documents within folder "notes_TwoUserCase"
    '''
    return rr.intersection([MAC_capacity_independentMessages(P1, P2, g13, g23),
                            MAC_capacity_independentMessages(P1, P2, g14, g24)
                            ])

def Point2Point_capacitiy(P, h=1):
    '''
    Capacity region of a point-to-point channel.
    Y = h*X + Z
    R <= 0.5*log2( 1+(h^2)*P ) . 
    '''      

def FCo_scheme(P3, P4):    
    '''
    Full-Cooperation 2nd hop, 
    two relays which have both been given (w1, w2), will both send both messages 
    using the coherent codewords:  
    Y = X3 + X4 + Z 
      = X3 + ( sqrt(P4/P3) )*X3 + Z  # full cooperation!
      = ( 1 + sqrt(P4/P3) ) * X3 + Z 
    '''
    rSum = _computeC(P3+P4+2*np.sqrt(P3*P4 ) )
    inequalities = [(1, 1, rSum) ] # R1 + R2 <= rSum 
    return  rr.InequalityRegion(inequalities)


def NoInterference_1stHop(P1, P2):
    '''
    the intersection of the capacity regions of two point-2-point channels
    '''
    r1 = _computeC( P1 )
    r2 = _computeC( P2 )
    inequalities = [(1, 0, r1), # R1 + 0  <= r1
                    (0, 1, r2)] #  0 + R2 <= r2  
    return  rr.InequalityRegion(inequalities)

def NoInterference_2ndHop(P3, P4):
    '''
    is actually the capacity region of a two-user MAC with independent messages
    '''
    return MAC_capacity_independentMessages(P3, P4)
    

def TwoHops_CF_ICFsch1(Ps, P3, P4, g13, g14, g23, g24):
    '''
    First, take the intersection of rate regions of 1st and 2nd hop
    Then, the the convex hull.
    Note: 
    Conceptually, we should do 
    1). for bigR1, take the intersection of CF_bigR1 and ICF_bigR2
    2). for bigR2, take the intersection of CF_bigR2 and ICF_bigR2
    3). take the union of the above two intersections 
    4). take the convex hull. 
    Actually, we did
    a) take the union of ICF_bigR1 and ICF_bigR2
    b) take the intersection of a) and CF_whole
    c) take the convex hull of b)
    '''
    return rr.intersection( [rr.union( [ICF_sch1_bigR1(P3, P4), 
                                 ICF_sch1_bigR2(P3, P4)]),
                      CF_scheme_integerCoeff(Ps, g13, g14, g23, g24)
                      ] ).convexHull()

def TwoHops_CF_ICFsch1_fromRegions(icf_sch1_bigR1, 
                                   icf_sch1_bigR2,
                                   cf_scheme_integerCoeff):
    '''
    compute the achievable rate region for the (whole) two-hop network from the 
    given rate regions for the 1st and 2nd hops 
    '''
    return rr.intersection( [rr.union( [icf_sch1_bigR1,
                                        icf_sch1_bigR2]),
                             cf_scheme_integerCoeff
                             ] ).convexHull()                      


def TwoHops_CF_ICFsch2(Ps, P3, P4, g13, g14, g23, g24):
    '''
    First, take the intersection of rate regions of 1st and 2nd hop
    Then, the the convex hull.
    Note: 
    Conceptually, we should do 
    1). for bigR1, take the intersection of CF_bigR1 and ICF_bigR2
    2). for bigR2, take the intersection of CF_bigR2 and ICF_bigR2
    3). take the union of the above two intersections 
    4). take the convex hull. 
    Actually, we did
    a) take the union of ICF_bigR1 and ICF_bigR2
    b) take the intersection of a) and CF_whole
    c) take the convex hull of b)
    '''
    return rr.intersection( [rr.union( [ICF_sch2_bigR1(P3, P4), 
                                 ICF_sch2_bigR2(P3, P4)]),
                      CF_scheme_integerCoeff(Ps, g13, g14, g23, g24)
                      ] ).convexHull()
                      
def TwoHops_CF_ICFsch2_fromRegions(icf_sch2_bigR1, 
                                   icf_sch2_bigR2,
                                   cf_scheme_integerCoeff):
    '''
    compute the achievable rate region for the (whole) two-hop network from the 
    given rate regions for the 1st and 2nd hops 
    '''
    return rr.intersection( [rr.union( [icf_sch2_bigR1,
                                        icf_sch2_bigR2]),
                             cf_scheme_integerCoeff ] ).convexHull()                        


def TwoHops_CF_ICFsch3(Ps, P3, P4, g13, g14, g23, g24):
    '''
    First, take the intersection of rate regions of 1st and 2nd hop
    Then, the the convex hull.
    Note: 
    Conceptually, we should do 
    1). for bigR1, take the intersection of CF_bigR1 and ICF_bigR2
    2). for bigR2, take the intersection of CF_bigR2 and ICF_bigR2
    3). take the union of the above two intersections 
    4). take the convex hull. 
    Actually, we did
    a) take the union of ICF_bigR1 and ICF_bigR2
    b) take the intersection of a) and CF_whole
    c) take the convex hull of b)
    '''
    alternative = True
    if alternative:
        # use alternative_rate_region.py 
        return rr.intersection( [rr.union( [alternative_compute.ICF_sch3_bigR1(P3, P4), 
                                 alternative_compute.ICF_sch3_bigR2(P3, P4)]),
                      CF_scheme_integerCoeff(Ps, g13, g14, g23, g24)] ).convexHull()
    else:
        # use rate_region.py
        return rr.intersection( [rr.union( [ICF_sch3_bigR1(P3, P4), 
                                 ICF_sch3_bigR2(P3, P4)]),
                      CF_scheme_integerCoeff(Ps, g13, g14, g23, g24) ] ).convexHull()


def TwoHops_CF_ICFsch3_fromRegions(cf_scheme_integerCoeff,
                                   *icf_sch3
                                   ):
    '''
    compute the achievable rate region for the (whole) two-hop network from the 
    given rate regions for the 1st and 2nd hops 
    '''
    if len(icf_sch3) == 1:
        return rr.intersection( [rr.union( [icf_sch3[0],
                                            icf_sch3[0].symmetricRegion_yEqualToX() ]),
                                 cf_scheme_integerCoeff] ).convexHull() 
    elif len(icf_sch3) == 2: 
        return rr.intersection( [rr.union( [icf_sch3[0],
                                            icf_sch3[1]]),
                                 cf_scheme_integerCoeff] ).convexHull()  

def TwoHops_DFIntegerCoeff_FCo(Ps, P3, P4, g13, g14, g23, g24):
    '''
    Both relays decode both messages and then forward these two messages 
    coherently to the destination node.
    '''
    return rr.intersection( [DF_scheme(Ps, Ps, g13, g14, g23, g24),
                      FCo_scheme(P3, P4) ]).convexHull()
                      

def TwoHops_DFIntegerCoeff_FCo_fromRegions(df_scheme, fCo_scheme):
    '''
    take the convex hull of the intersection of two rate regions for 1st and 
    2nd hop. 
    '''
    return rr.intersection( [df_scheme, 
                             fCo_scheme ]).convexHull()

def TwoHops_noInterference(P1, P2, P3, P4):
    '''
    take the convex hull of the intersection. 
    '''
    return rr.intersection( [NoInterference_1stHop(P1, P2),
                             NoInterference_2ndHop(P3, P4) ]).convexHull()


def TwoHops_noInterference_fromRegions(noInterference_1stHop, noInterference_2ndHop):
    '''
    take the convex hull of the intersection. 
    '''
    return rr.intersection( [noInterference_1stHop,
                             noInterference_2ndHop  ] ).convexHull()
                             











    










