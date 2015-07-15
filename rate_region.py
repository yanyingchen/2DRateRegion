'''
Created on Mar 27, 2015

Compute the rate region, specified by a collection of linear inequalities 
x*c0+y*c1<=c2
in the first quadrant. 

Note that k=50000, i.e. the square (0,0) to (k,k), is used as bounding box. 
I.e., we assume that no rate is larger than k. 
@author: yanying
'''

import shapely.geometry as geo

# for plotting
import numpy as np
import matplotlib.pyplot as pyplot

class RateRegion(object):
    '''
    A 2D set of points in the plane
    '''        
    # the size of boundary so that we have finite regions, which are used for 
    # denoting an half plane
    k = 50000          
    
    def __init__(self,geometry):
        '''
        Constructs a region given its geometry
        '''
        self._geometry = geometry
    
    def convexHull(self):
        '''
        Returns the convex hull of this region
        '''
        return RateRegion(self._geometry.convex_hull)
    
    def symmetricRegion_yEqualToX(self):
        '''
        Returns the symmetric region w.r.t. line y=x of this region
        '''
        boundaryPoints = np.asarray(self._geometry.boundary)
        # switch the (x,y) coordinates
        symmetricBoundary = np.vstack((boundaryPoints[:,1], 
                                       boundaryPoints[:, 0] )).transpose()
        return RateRegion( geo.MultiPoint(symmetricBoundary).convex_hull)                                        
    
    def _repr_svg_(self):
        '''
        Interface for rich displays in ipython
        '''
        return self.geometry()._repr_svg_()
    
    def geometry(self):
        firstQuadrantBound = geo.box(0, 0, RateRegion.k, RateRegion.k)
        self._geometry = self._geometry.intersection(firstQuadrantBound)
        return self._geometry
    
def _setOp(op,regionList):
    '''
    Returns a region that is the combination of the regions given in a list
    The actual combination is given by op
    '''
    geometry = regionList[0]._geometry
    for r in regionList[1:]:
        gop = op.__get__(geometry)
        geometry = gop(r._geometry)
    return RateRegion(geometry)    

def unionAndConvexhull(regionList):
    '''
    Returns the convex hull of the union of the regions given in a list.
    In the union step, we speed up by just keeping all boundary-points of the
    given regions, rather than applying union operation directly on regions for 
    many times
    '''
    points = [np.asarray(r._geometry.boundary ) for r in regionList]
    geometry = geo.MultiPoint( np.vstack(points ) ).convex_hull
    return RateRegion( geometry )
        
      
     
def union(regionList):
    '''
    Returns a region that is the union of the regions given in a list
    '''
    return _setOp(geo.base.BaseGeometry.union,regionList)

def intersection(regionList):
    '''
    Returns a region that is the intersection of the regions given in a list
    '''
    return _setOp(geo.base.BaseGeometry.intersection,regionList)


def plot(rr,*args,**kwargs):
    '''
    Plot a region in the current figure of matplotlib
    '''
    u = np.asarray(rr.geometry().boundary)
    if 'axes' in kwargs:
        return kwargs['axes'].plot(u[:,0],u[:,1],*args,**kwargs)
    else:
        return pyplot.plot(u[:,0],u[:,1],*args,**kwargs)
        
class InequalityRegion(RateRegion):
    '''
    A region consisting on the intersection of multiple regions, typically,
    HalfPlane's
    '''
    
    def _oneInequalityRegion(self,c0,c1,c2):
        '''
        Compute the rate region for one inequality: x*c0+y*c1<=c2
        
        We achieve this by first building a geometry that extents outside of
        the boundary and then trimming it (via intersection) to the
        target square
        '''
        
        k = 2*RateRegion.k # twice the size of boundary
    
        if abs(c1) > abs(c0):
            # separating lines that are "mostly horizontal" abs(c1)>abs(c0)
            # this guarantees in particular abs(c1) > 0
            p0 = (-k , (c2+k*c0)/c1 )
            p1 = ( k , (c2-k*c0)/c1 ) 
            if c1>0: # valid region is below the separating line
                p2 = (-k,-k)
                p3 = ( k,-k) 
            else: # valid region is above the separating line
                p2 = (-k, k)
                p3 = ( k, k) 
        else:
            # separating lines that are "mostly vertical" abs(c1)<abs(c0)
            # this guarantees in particular abs(c0) > 0
            p0 = ( (c2+k*c1)/c0 , -k )
            p1 = ( (c2-k*c1)/c0 ,  k )
            if c0>0: # valid region is to the left the separating line
                p2 = (-k,-k)
                p3 = (-k, k) 
            else: # valid region is to the right the separating line
                p2 = (k,-k)
                p3 = (k, k)   
                
        geometry = geo.MultiPoint([p0,p1,p2,p3]).convex_hull #geometry extends outside the target square
        return geometry   
    
    def __init__(self, params):
        '''
        Construct a RateRegion which is the intersection of several half planes
        params is the list of planes to intersect, where each plan is specified 
        in the canonical form [c0, c1, c2]: x*c0+y*c1<=c2 . 
        '''
        
        self._geometry = self._oneInequalityRegion(*params[0])
        for r in params[1:]:
            self._geometry = self._geometry.intersection(self._oneInequalityRegion(*r))

            