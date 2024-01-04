# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 21:10:21 2022

@author: HP
"""

"""
Script that contains the function required for lens optimisation.
"""

# MODULES
import scipy.optimize as op
import ray_class as rc
import optical_element_class as oec
import matplotlib.pyplot as plt
from copy import deepcopy

#FUNCTIONS
def lens_optimise(rays, intercept, seperation, n, aperture_radius, 
                  focus):
    """
    Function that uses scipy.minimise to vary the curvatures of a two surface 
    lens such that the RMS spread at the output plane is minimised.
    
    Returns the minimised RMS spread after fun: and the curvatures of the 
    surfaces are given as an array after x: . The first value in the array 
    corresponds to the curvature of the leftmost surface.
    
    WARNING: Function can take a long time to run with a large number of rays.
    
    Parameters
    ----------
    rays : list
        List of rays travelling parallel to the optical axis.
    intercept : float
        The z-coordinate of the intercept of the leftmost surcace with the z-
        axis.
    seperation : float
        Seperation between the two lenses of the system.
    n : float
        Refractive index of the lens (function assumes surround refractive
        index = 1).
    aperture_radius : float
        Maximum perpendicular distance from the optical axis to the top of the
        aperture.
    focus : array_like
        List that contains the lower bound and upper bound for the range in
        which the lens aims to focus the rays to.
    
    """
    rays1 = rays
    # Mean focal distance between the lens and the max/min acceptable focal
    # positions.
    focus_mean = (focus[0]+focus[1])/2
    focal_distance = focus_mean - intercept
    
    # Theoretical equation for the radius of a plano convex lens given the 
    # focal distance.
    radius = (focal_distance)*(n-1)
    curvature = 1/radius
    
    # Using the curvatures of a plano-convex setup lens as an initial guess 
    # for the optimise function.
    initial_guess = [curvature,0]
    
    def best_form(curves):
        """
        Function that takes a list of two curvatures as its parameter and
        returns RMS radius. This is the function that is minimised by the
        scipy minimize function.
        """
        # Position of the z-intecept of the left lens.
        curvature_1 = curves[0]
        curvature_2 = curves[1]
        Z = intercept
        ray_test= rc.Ray([0,0.1,0],[0,0,1])
        s1 = oec.SphericalRefraction(Z,curvature_1,1,n,aperture_radius)
        s2 = oec.SphericalRefraction(Z+seperation,curvature_2,n,1,
                                     aperture_radius)
        
        # Propogating paraxial ray through the two surfaces to find the z 
        # position of the actual paraxial focus of the two setups.
        s1.propogate_ray(ray_test)
        s2.propogate_ray(ray_test)
        focal_point = ray_test.focal_z()
        out = oec.OutputPlane(focal_point)
        
        # Deepcopy used so that the original rays are used for each iteration
        # of the minimise function.
        rays = deepcopy(rays1)
        rms = oec.ray_trace_2d(rays,[s1,s2],out,rms=True,spot=True)
        plt.show()
        print("Focal point: %s mm "%(focal_point))
        
        return rms
    
    def constraint(curves):
        """
        Function that is used as a constraint such that the ray is 
        focussed to a specific point when optimising.
        """
        curvature_1 = curves[0]
        curvature_2 = curves[1]
        
        # Using a paraxial test ray to find the current paraxial focal 
        # position.
        ray_test= rc.Ray([0,0.0001,0],[0,0,1])
        s1 = oec.SphericalRefraction(intercept,curvature_1,1,n, 
                                     aperture_radius)
        s2 = oec.SphericalRefraction(intercept+seperation,curvature_2,n,1,
                                     aperture_radius)
        
        # Propogating paraxial ray through the two surfaces to find the z 
        # position of the paraxial focus of the two setups.
        s1.propogate_ray(ray_test)
        s2.propogate_ray(ray_test)
        focal_point = ray_test.focal_z()
        
        # Focal point should = designated focus as a constraint.
        return focal_point
    
    
    def constraint1(curves):
        """Constraint for the minimum focal position boundary."""
        delta = constraint(curves) - focus[0]
        return delta
   
    def constraint2(curves):
        """Constraing for the maximum focal position boudnary."""
        delta = -(constraint(curves) - focus[1])
        return delta
    
    con1 = {'type':'ineq','fun':constraint1}
    con2 = {'type':'ineq','fun':constraint2}
    cons = [con1,con2]
    rms_optimise = op.minimize(best_form,initial_guess, constraints = cons)
    
    return rms_optimise
