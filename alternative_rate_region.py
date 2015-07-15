'''
Created on Apr 5, 2015

Please first check rate_region.py 
This alternative file is developed to speed up the computation.
*******************************
**Created on Mar 27, 2015

Note that as far as the simulations in ICFjournal is concerned, the inequalities 
that define a rate region can only of the following formats: 
    (a)  1*x + 0*y <= xmax   
    (b)  0*x + 1*y <= ymax   
    (c)  1*x + 1*y <= xysum  
    (d1) 1*x + (-1)*y <= 0   
    (d2) (-1)*x + 1*y <= 0   
Thus, the implementation here applies only to these inequality types. 
*******************************
@author: yanying
'''





import shapely.geometry as geo

from scipy.spatial import ConvexHull

# for plotting
import numpy as np
import matplotlib.pyplot as plt
 
# to make make the first quadrant finite, only consider the regions inside 
# square (0,0) - (k,k) 
k = 100000       



def _compute_yes_yes_yes_noo_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (a)  1*x + 0*y <= xmax   
        (b)  0*x + 1*y <= ymax   
        (c)  1*x + 1*y <= xysum  
    ''' 
    if (xmax<xysum):
        if (ymax<xysum):
            # the region is due to (a, b, c)
            if (xmax+ymax) < xysum:
                # the region is a rectangle
                return np.array([[0,0], [xmax,0], [xmax,ymax], [0,ymax] ] )
            else:
                # the region is a pentagon 
                return np.array([[0,0], [xmax,0], [xmax,xysum-xmax], 
                         [xysum-ymax,ymax], [0,ymax]   ])     
        else:
            # the region is due to (a, c)
            return np.array([[0,0], [xmax,0], [xmax,xysum-xmax], [0,xysum] ])
    else:
        # the region is due to (b) and/or (c)
        if (ymax<xysum):
            # the region is due to (b, c)
            return  np.array([[0,0], [xysum,0], [xysum-ymax, ymax], [0, ymax] ]) 
        else:
            # the region is due to (c)
            return np.array([[0,0], [xysum,0], [0,xysum] ] )


def _compute_yes_yes_noo_noo_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (a)  1*x + 0*y <= xmax   
        (b)  0*x + 1*y <= ymax 
    ''' 
    # the region is a rectangle
    return np.array([[0,0], [xmax,0], [xmax,ymax], [0,ymax] ])    


def _compute_yes_noo_yes_noo_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (a)  1*x + 0*y <= xmax   
        (c)  1*x + 1*y <= xysum  
    ''' 
    if (xmax<xysum):
        # the region is due to (a, c)
        return np.array([[0,0], [xmax,0], [xmax,xysum-xmax], [0,xysum] ])
    else:
        # the region is due to (c)
        return np.array([[0,0], [xysum,0], [0,xysum] ] )
    
    

def _compute_yes_noo_noo_noo_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (a)  1*x + 0*y <= xmax   
    ''' 
    # square (0,0)-(k,k) is the bounding area
    return np.array([[0,0], [xmax,0], [xmax,k], [0,k] ])


def _compute_noo_yes_yes_noo_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (b)  0*x + 1*y <= ymax   
        (c)  1*x + 1*y <= xysum    
    ''' 
    if (ymax<xysum):
        # the region is due to (b, c)
        return  np.array([[0,0], [xysum,0], [xysum-ymax, ymax], [0, ymax] ]) 
    else:
        # the region is due to (c)
        return np.array([[0,0], [xysum,0], [0,xysum] ] )


def _compute_noo_yes_noo_noo_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (b)  0*x + 1*y <= ymax   
    '''
    # square (0,0)-(k,k) is the bounding area
    return np.array([[0,0], [k,0], [k,ymax], [0,ymax] ]) 
     

def _compute_noo_noo_yes_noo_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (c)  1*x + 1*y <= xysum 
    '''
    return np.array([[0,0], [xysum,0], [0,xysum] ] )


def _compute_noo_noo_noo_noo_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        - - no constraints - - 
    '''
    # square (0,0)-(k,k) is the bounding area
    return np.array([[0,0], [k,0], [k,k], [0,k] ] )


def _compute_yes_yes_yes_noo_yes(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (a)  1*x + 0*y <= xmax   
        (b)  0*x + 1*y <= ymax   
        (c)  1*x + 1*y <= xysum  
        (d2) 1*x + (-1)*y <= 0  when yBigger is True
    '''
    # consider the region due to (c, d2) as the base region. 
    # The following discussion is on (a,b).
    if (xmax>=(0.5*xysum)) and (ymax>=xysum):
        # remove constraints (a) and (b)
        return np.array([[0,0], [0.5*xysum, 0.5*xysum], [0, xysum] ])
    if (xmax>=(0.5*xysum)) and (ymax<xysum):
        # remove constraint (a). Only need to consider (b, c, d2)
        if ymax>(0.5*xysum):
            return np.array([[0,0], [0.5*xysum,0.5*xysum], [xysum-ymax,ymax], [0,ymax] ])
        else:
            return np.array([[0,0], [ymax,ymax], [0,ymax] ])
    if (xmax<(0.5*xysum)) and (ymax>=xysum):
        # consider constraints (a, c, d2)
        return np.array([[0,0], [xmax, xmax], [xmax, xysum-xmax], [0,xysum] ])
    if (xmax<(0.5*xysum)) and (ymax<xysum):
        # consider the intersection points of 
        # (a and c): (xmax, xysmum-xmax)
        # (a and d2): (xmax, xmax)
        if ymax>(xysum-xmax):
            return np.array([[0,0], [xmax, xmax], [xmax, xysum-xmax],
                             [xysum-ymax, ymax], [0,ymax] ])
        elif ymax==(xysum-xmax):
            return np.array([[0,0], [xmax, xmax], [xmax, ymax], [0,ymax] ])
        elif ymax>xmax:
            return np.array([[0,0], [xmax,xmax], [xmax,ymax], [0,ymax]  ])
        elif ymax==xmax:
            return np.array([[0,0], [ymax,ymax], [0,ymax] ])
        else: #ymax<xmax
            return np.array([[0,0], [ymax,ymax], [0,ymax]  ])
    

def _compute_yes_yes_noo_noo_yes(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (a)  1*x + 0*y <= xmax   
        (b)  0*x + 1*y <= ymax   
        (d2) 1*x + (-1)*y <= 0  when yBigger is True
    '''
    # consider the triangle defined by (a, d2), denoted as base region
    # discuss the relative position of line x=xmax and line x=ymax
    if xmax>=ymax:
        # return the base region
        return np.array([[0,0], [ymax,ymax], [0,ymax] ])
    else:
        return np.array([[0,0], [xmax,xmax], [xmax,ymax], [0,ymax] ])
    

def _compute_yes_noo_yes_noo_yes(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (a)  1*x + 0*y <= xmax   
        (c)  1*x + 1*y <= xysum  
        (d2) 1*x + (-1)*y <= 0  when yBigger is True
    '''
    # consider the region due to (c, d2) as the base region.
    if xmax<(0.5*xysum):
        return np.array([[0,0], [xmax,xmax], [xmax,xysum-xmax], [0,xysum] ])
    else:
        return np.array([[0,0], [0.5*xysum, 0.5*xysum], [0,xysum] ])


def _compute_yes_noo_noo_noo_yes(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (a)  1*x + 0*y <= xmax   
        (d2) 1*x + (-1)*y <= 0  when yBigger is True
    '''
    return np.array([[0,0], [xmax,0], [xmax,xmax] ])


def _compute_noo_yes_yes_noo_yes(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (b)  0*x + 1*y <= ymax   
        (c)  1*x + 1*y <= xysum  
        (d2) 1*x + (-1)*y <= 0  when yBigger is True
    '''
    # consider the rgion due to (c, d2) as the base region
    # horizontal lines: y=xysum, 0.5*xysum
    if ymax>=xysum:
        return np.array([[0,0], [0.5*xysum,0.5*xysum], [0,xysum] ])
    elif ymax>(0.5*xysum):
        return np.array([[0,0], [xysum-ymax, xysum-ymax], 
                         [xysum-ymax,ymax], [0,ymax] ])
    else:
        return np.array([[0,0], [ymax,ymax], [0,ymax] ])
        

def _compute_noo_yes_noo_noo_yes(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (b)  0*x + 1*y <= ymax   
        (d2) 1*x + (-1)*y <= 0  when yBigger is True
    '''
    return np.array([[0,0], [ymax,ymax], [0,ymax] ])


def _compute_noo_noo_yes_noo_yes(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (c)  1*x + 1*y <= xysum     
        (d2) 1*x + (-1)*y <= 0  when yBigger is True
    '''
    return np.array([[0,0], [0.5*xysum, 0.5*xysum], [0,xysum] ])


def _compute_noo_noo_noo_noo_yes(xmax, ymax, xysum, xBigger, yBigger):
    '''
    Returns a numpy array that matches the boundary for the region given by 
        (d2) 1*x + (-1)*y <= 0  when yBigger is True
    '''
    # square (0,0)-(k,k) is the bounding area
    return np.array([[0,0], [k,k], [0,k] ])


def _compute_yes_yes_yes_yes_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    '''
    # swap the x-y coordinate system
    swappedArray = _compute_yes_yes_yes_noo_yes(ymax, xmax, xysum, yBigger, xBigger)
    # swap (x,y) again to go back 
    return swappedArray[:, [1,0]]


def _compute_yes_yes_noo_yes_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    '''
    # swap the x-y coordinate system
    swappedArray = _compute_yes_yes_noo_noo_yes(ymax, xmax, xysum, yBigger, xBigger)
    # swap (x,y) again to go back 
    return swappedArray[:, [1,0]]


def _compute_yes_noo_yes_yes_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    '''
    # swap the x-y coordinate system
    swappedArray = _compute_yes_noo_yes_noo_yes(ymax, xmax, xysum, yBigger, xBigger)
    # swap (x,y) again to go back 
    return swappedArray[:, [1,0]]


def _compute_yes_noo_noo_yes_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    '''
    # swap the x-y coordinate system
    swappedArray = _compute_yes_noo_noo_noo_yes(ymax, xmax, xysum, yBigger, xBigger)
    # swap (x,y) again to go back 
    return swappedArray[:, [1,0]]


def _compute_noo_yes_yes_yes_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    '''
    # swap the x-y coordinate system
    swappedArray = _compute_noo_yes_yes_noo_yes(ymax, xmax, xysum, yBigger, xBigger)
    # swap (x,y) again to go back 
    return swappedArray[:, [1,0]]


def _compute_noo_yes_noo_yes_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    '''
    # swap the x-y coordinate system
    swappedArray = _compute_noo_yes_noo_noo_yes(ymax, xmax, xysum, yBigger, xBigger)
    # swap (x,y) again to go back 
    return swappedArray[:, [1,0]]


def _compute_noo_noo_yes_yes_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    '''
    # swap the x-y coordinate system
    swappedArray = _compute_noo_noo_yes_noo_yes(ymax, xmax, xysum, yBigger, xBigger)
    # swap (x,y) again to go back 
    return swappedArray[:, [1,0]]


def _compute_noo_noo_noo_yes_noo(xmax, ymax, xysum, xBigger, yBigger):
    '''
    '''
    # swap the x-y coordinate system
    swappedArray = _compute_noo_noo_noo_noo_yes(ymax, xmax, xysum, yBigger, xBigger)
    # swap (x,y) again to go back 
    return swappedArray[:, [1,0]]
'''
Note that
    [xmax, ymax, xysum, yBigger], where 
    1: xmax, ymax, xysum are either positive numbers or Inf
    2: yBigger is a boolean variable         
    3: The value of [xmax, ymax, xysum, yBigger] indicates a set of inequalities, 
    which contains one or several of the following inequalities   
    (a)  1*x + 0*y <= xmax   
    (b)  0*x + 1*y <= ymax   
    (c)  1*x + 1*y <= xysum  
    (d1) (-1)*x + 1*y <= 0  when xBigger is True
    (d2) 1*x + (-1)*y <= 0  when yBigger is True
'''
# put to following dictionary inside the class?? outside of the class? 
yes = True # yes, impose the inequality constraint
noo = False # no, do not impose the inequality constraint
_toNumpyArrayDict = {(yes, yes, yes, noo, noo): _compute_yes_yes_yes_noo_noo,
                     (yes, yes, noo, noo, noo): _compute_yes_yes_noo_noo_noo,
                     (yes, noo, yes, noo, noo): _compute_yes_noo_yes_noo_noo,
                     (yes, noo, noo, noo, noo): _compute_yes_noo_noo_noo_noo,
                     (noo, yes, yes, noo, noo): _compute_noo_yes_yes_noo_noo,
                     (noo, yes, noo, noo, noo): _compute_noo_yes_noo_noo_noo,
                     (noo, noo, yes, noo, noo): _compute_noo_noo_yes_noo_noo,
                     (noo, noo, noo, noo, noo): _compute_noo_noo_noo_noo_noo,
                     
                     (yes, yes, yes, noo, yes): _compute_yes_yes_yes_noo_yes,
                     (yes, yes, noo, noo, yes): _compute_yes_yes_noo_noo_yes,
                     (yes, noo, yes, noo, yes): _compute_yes_noo_yes_noo_yes,
                     (yes, noo, noo, noo, yes): _compute_yes_noo_noo_noo_yes,
                     (noo, yes, yes, noo, yes): _compute_noo_yes_yes_noo_yes,
                     (noo, yes, noo, noo, yes): _compute_noo_yes_noo_noo_yes,
                     (noo, noo, yes, noo, yes): _compute_noo_noo_yes_noo_yes,
                     (noo, noo, noo, noo, yes): _compute_noo_noo_noo_noo_yes,
                    
                     (yes, yes, yes, yes, noo): _compute_yes_yes_yes_yes_noo,
                     (yes, yes, noo, yes, noo): _compute_yes_yes_noo_yes_noo,
                     (yes, noo, yes, yes, noo): _compute_yes_noo_yes_yes_noo,
                     (yes, noo, noo, yes, noo): _compute_yes_noo_noo_yes_noo,
                     (noo, yes, yes, yes, noo): _compute_noo_yes_yes_yes_noo,
                     (noo, yes, noo, yes, noo): _compute_noo_yes_noo_yes_noo,
                     (noo, noo, yes, yes, noo): _compute_noo_noo_yes_yes_noo, 
                     (noo, noo, noo, yes, noo): _compute_noo_noo_noo_yes_noo
}  

class ICFsch3RateRegion(object):
    '''
    A region consisting on the intersection of multiple regions, typically,
    of type (a), (b), (c), (d1) or (d2) 
    The internal representation of  a universal format [xmax, ymax, xysum, yBigger]
    '''
    

    
    def _oneUniversalRegion(self,c0,c1,c2):
        '''
        Compute the rate region for one inequality: x*c0+y*c1<=c2
                
        We save this inequality region in our "universal form",represented by 
        [xmax, ymax, xysum, xBigger, yBigger], where 
        1: xmax, ymax, xysum are either positive numbers or Inf
        2: xBigger, yBigger are a boolean variable
        3. xBigger and yBigger cannot be true at the same time (usually)         
        3: The value of [xmax, ymax, xysum, xBigger, yBigger] indicates a set of
        inequalities, which contains one or several of the following inequalities   
        (a)  1*x + 0*y <= xmax   
        (b)  0*x + 1*y <= ymax   
        (c)  1*x + 1*y <= xysum  
        (d1) (-1)*x + 1*y <= 0  when xBigger is True
        (d2) 1*x + (-1)*y <= 0  when yBigger is True
        '''
        case = c0 + c1 
        if case==2:
            # 1*x + 1*y <= c2
            xmax, ymax, xysum, xBigger, yBigger = float('Inf'), float('Inf'), c2, False, False
        elif case==1: 
            if c0==1:
                # 1*x + 0*y <= c2
                xmax, ymax, xysum, xBigger, yBigger = c2, float('Inf'), float('Inf'), False, False
            else:
                # 0*x + 1*y <= c2
                xmax, ymax, xysum, xBigger, yBigger = float('Inf'), c2, float('Inf'), False, False
        else:
            if c0==1: 
                # 1*x + (-1)*y <= 0
                xmax, ymax, xysum, xBigger, yBigger = float('Inf'), float('Inf'), float('Inf'), False, True
            else:
                # (-1)*x + 1*y <== 0
                xmax, ymax, xysum, xBigger, yBigger = float('Inf'), float('Inf'), float('Inf'), True, False
        return [xmax, ymax, xysum, xBigger, yBigger]
    
   
     

    
    def __init__(self, 
                 listOfCanonicalInequalities = None,
                 listOfUniversalInequalities = None):
        '''
        Construct a RateRegion which is the intersection of several inequalities.
        - listOfCanonicalInequalities is the list of inequalities to intersect, where each inequality 
        is specified in the canonical form [c0, c1, c2]: x*c0+y*c1<=c2 . 
        - listOfUniversalInequalities is the list inequalities in their universal representation: 
        [xmax, ymax, xysum, yBigger]
         
        '''
        if listOfCanonicalInequalities and (not listOfUniversalInequalities) :
            listOfUniversalInequalities = [self._oneUniversalRegion(c0, c1, c2) for (c0, c1, c2) in listOfCanonicalInequalities]
        elif listOfUniversalInequalities and (not listOfCanonicalInequalities):
            pass            
        else:
            raise Exception('please use only one parameter: \n - use listOfCanonicalInequalities with or without the keyword \n or \n - use listOfUniversalInequalities WITH the keyword')
        
        self._inequalities = ICFsch3RateRegion.intersection(listOfUniversalInequalities) 
        # self.computeBoundaryPoints()
            
    def computeBoundaryPoints(self):
        '''
        compute _asNumpyArray variable that represents the polygon around this region. 
        
        Note that the points are listed counter clockwise and the first one is at
        the bottom left corner.  
        '''
        xmax, ymax, xysum, xBigger, yBigger = self._inequalities
        self._asNumpyArray = _toNumpyArrayDict[(xmax!=float('Inf'), 
                                 ymax!=float('Inf'), 
                                 xysum!=float('Inf'),
                                 xBigger, 
                                 yBigger)](xmax, ymax, xysum, xBigger, yBigger)
        
    
    @staticmethod
    def intersection(listOfUniversalInequalities):
        '''
        Returns the intersection of a list of universal regions
        
        Note that one universal region is [xmax, ymax, xysum, xBigger, yBigger]
        '''
        asArray = np.asarray(listOfUniversalInequalities)        
        # The inequality with smallest right hand side value dominates.
        # take the minimum among (real numbers, Inf). 
        xmax, ymax, xysum = np.amin(asArray[:, 0:3], axis=0)
        # Impose the constraint as long as there exists at least one constraint  
        # take the maximum among (True, False) 
        xBigger, yBigger = np.amax(asArray[:, 3:5], axis=0)          
        return [xmax, ymax, xysum, xBigger, yBigger]
                                   
def intersection(regionList):
    '''
    Returns a region that is the intersection of the regions given in a list
    '''
    return ICFsch3RateRegion(listOfUniversalInequalities = [r._inequalities for r in regionList]) 

def unionAndConvexhull(regionList):
    '''
    Returns a numpy array, which is the convex hull of the union of the regions 
    given in a list as a numpy array.
    In the union step, we speed up by just keeping all boundary-points of the
    given regions, rather than applying union operation directly on regions for 
    many times
    '''
    for r in regionList:
        r.computeBoundaryPoints()
    points = np.vstack([r._asNumpyArray for r in regionList])
    hull = ConvexHull(points)
    # oneRateRegion = RateRegion(listOf)
    return points[hull.vertices,:]
    
            