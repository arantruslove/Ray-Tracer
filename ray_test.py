# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 17:31:35 2022

@author: HP
"""
"""Script used to test the code produced for each task."""

# No need for testing task 1.

# %% MODULES
import optical_element_class as oec
import ray_class as rc
import lens_optimisation as lop
import matplotlib.pyplot as plt

# %% Task 2
"""
Task 2: Testing the behaviour of the ray class by testing initialisation and
appending points to the ray.
"""
# Creating a random instance of the Ray class.
r2 = rc.Ray(p=[0, 0, 2], k=[1, 2, 1])
print("Vertices of the random ray before appending:%s" % (r2.vertices()))

# Testing the append method of the Ray class.
r2.append(p=[1, 1, 1], k=[0, 0, 1])
print("Vertices of the random ray after appending:%s" % (r2.vertices()))

# Testing the initialisation of the class.
r2.__init__(p=[1, 1, 1], k=[0, 0, 1])
print("Vertices of the random ray after initialising:%s" % (r2.vertices()))

# %% Task 3
""" 
Task 3: Testing the behaviour of the Spherical_Refraction class which is an
    instance of the Optical_Element class.
"""
# Instance of the SphericalRefraction class.
s3 = oec.SphericalRefraction(3, 0.1, 1, 1.2, 30)
print(
    "Coordinates of the centre of the spherical refracting surface:%s" % (s3.centre())
)

# Initialising the class with a negative curvature.
s3.__init__(3, -1 / 10, 1, 1.2, 30)
print(
    """Coordinates of the centre of the spherical refracting surface with a 
negative curvature:%s"""
    % (s3.centre())
)

# %% Task 4
"""
Task 4, part 1: Testing an intercept method of the Spherical_Refraction class 
which calculates the coordinates of the first valid intercept of the ray with 
the spherical suface.
"""
r4 = rc.Ray(p=[0, 0, 1], k=[1, 1, 1])
s4 = oec.SphericalRefraction(3, 0.1, 1, 1.2, 100)

# Finding the first intercept of the ray with the surface.
print("Coordinates of intercept are %s" % (s4.intercept(r4)))

# Finding the intercept of the ray with a surface with a negative curvature.
s4 = oec.SphericalRefraction(3, -1 / 10, 1, 1.2, 10)
print("Coordinates of intercept with a negative curvature are %s" % (s4.intercept(r4)))

# Initialising the ray so that it travels vertically upwards and seeing if the
# intercept function returns None due to no intercept of ray with the
# surface.
r4.__init__(p=[1, 1, 1], k=[0, 1, 0])
s4.__init__(3, 1 / 10, 1, 1.2, 10)
print(s4.intercept(r4))

# %% Task 4
"""
Part 2: Testing continued...
"""
# Testing a ray that crosses with the 'imaginary' sphere above the aperture
# radius but then propogates through the 'real' lens.
r4.__init__(p=[0, 10, 3], k=[0, -0.5, 1])
s4.__init__(3, 1 / 10, 1, 1.2, 3)
print("Coordinates of the intercept are:%s" % (s4.intercept(r4)))

# Initialising the spherical surface with an aperture radius less than the
# height of the height of the ray and seeing whether the intercept function as
# returns None the ray travels over the aperture.
r4.__init__(p=[0, 5, 0], k=[0, 0, 1])
s4.__init__(0, 1 / 10, 1, 1.2, 4.99)
print(s4.intercept(r4))

# %% Task 4
"""
Part 3: Testing the intercept of a ray with a flat surface (C=0).
"""
s4 = oec.SphericalRefraction(3, 0, 1, 1, 100)
r4 = rc.Ray(p=[0, 10, 0], k=[0, 0, 1])
print("Coordinates of intercept of a ray with a flat surface:%s" % (s4.intercept(r4)))

# Testing for the case when the theoretical intercept point is greater than
# the aperture radius.
s4 = oec.SphericalRefraction(3, 0, 1, 1, 2)
print("Coordinates of intercept are:%s" % (s4.intercept(r4)))
print(s4.refracted_direction(r4))

# %% Task 5, Part 1
"""
Task 5: Testing a method of the SphericalRefraction class that implements
Snell's law.
"""
# Propogating a flat surface with n1 < n2 (C=0).
r5 = rc.Ray(p=[0, 0, 0], k=[0, 0.8, 1])
s5 = oec.SphericalRefraction(3, 0, 1, 2, 100)
print("Refracted direction for converging lens is:%s" % (s5.refracted_direction(r5)))

# Testing the TIR conditions with n2/n1 = 0.833.
s5.__init__(3, 0, 1.2, 1, 100)

# Approaching with an angle of 0.675 rad => sin(0.675) = 0.625
# Sin(theta) < n2/n1 so no TIR is expected to occur.
print("Refracted direction for diverging lens is:%s" % (s5.refracted_direction(r5)))

# Approaching with an angle of sin(theta) = 0.894 > 0.833 so TIR should occur.
r5.__init__(p=[0, 0, 0], k=[0, 2, 1])
print("Refracted direction for diverging lens is:%s" % (s5.refracted_direction(r5)))

# %% Task 5, Part 2
"""
Task 5: Testing continued...
"""
# Testing to check that a ray approaching with a large angle of incidence is
# not affected by the TIR condition when n1 < n2.
s5 = oec.SphericalRefraction(3, 0, 1, 1.2, 16)
r5 = rc.Ray(p=[0, 0, 0], k=[0, 5, 1])
print("Refracted direction is:%s" % (s5.refracted_direction(r5)))

# %% Task 6 and task 7
"""
Task 6 and task 7: Testing the SphericalRefraction method 'propogate_ray'.
"""
# Propogating a ray through multiple surfaces and then checking the vertices
# of the ray after each propogation.

s61 = oec.SphericalRefraction(3, 1 / 100, 1, 1.2, 100)
r61 = rc.Ray(p=[0, 0, 0], k=[0, 1, 1])
print("Vertices before propogating through a surface:%s" % (r61.vertices()))

# Represents a ray passing above a refracting surface (no intercept with it.)
s62 = oec.SphericalRefraction(8, 0, 1, 1.2, 3)
s62.propogate_ray(r61)

# Propogating through first surface.
s63 = oec.SphericalRefraction(15, 1 / 100, 1.2, 1, 100)
s63.propogate_ray(r61)
print("Vertices after propogating through first surface:%s" % (r61.vertices()))

# Propogating through second, flat, surface with n1 > n2.
s62 = oec.SphericalRefraction(3, 0, 1, 1, 100)

s62.openLogFile()

s62.propogate_ray(r61)
print("Vertices after propogating through a second, flat surface:%s" % (r61.vertices()))

# %% Task 8 and task 9
"""
Task 8 and task 9: Testing the plotting function "ray_trace_2d" which takes
rays, refracting surfaces and an output plane as its paramters.

Testing the plotting function for a specific instance of z0 = 100 mm,
curvature = 0.03 mm^-1, n1 = 1.0, n2 = 1.5 and output plane at z = 250 mm.
"""
s81 = oec.SphericalRefraction(100, 0.03, 1, 1.5, 100)
r81 = rc.Ray([0, 5, 0], [0, 0, 1])
out = oec.OutputPlane(250)

# Plotting the ray/lens system.
oec.ray_trace_2d([r81], [s81], out)

# %% Task 10, Part 1
"""
Task 10: Testing the instance of a ray travelling parallel and close to the 
optical axis (0.1 mm) and comparing this to the expected paraxial focus point
of a spherical surface.
"""
r101 = rc.Ray([0, 0.1, 0], [0, 0, 1])
s101 = s81 = oec.SphericalRefraction(100, 0.03, 1, 1.5, 0.2)
out101 = oec.OutputPlane(250)
oec.ray_trace_2d([r101], [s81], out101)
print("Ray crosses optical axis at z = 200 mm")

# %% Task 10, Part 2
"""
Task 10: Testing the interaction of multiple rays with a surface and testing
the 'ray_bundle' function that creates a list of rays.
"""
r102 = rc.Ray([0, 0, 0], [0, 0, 1])
r103 = rc.Ray([0, 1, 0], [0, 0, 1])
r104 = rc.Ray([0, 2, 0], [0, 0, 1])
r105 = rc.Ray([0, -1, 0], [0, 0, 1])
r106 = rc.Ray([0, -2, 0], [0, 0, 1])
s102 = oec.SphericalRefraction(100, 0.03, 1, 1.5, 3)
oec.ray_trace_2d([r102, r103, r104, r105, r106], [s102], out)
print()

# Testing a "ray_bundle" function that makes a list of rays with the same
# direction.
print("Now testing the same rays but with the 'ray_bundle' function.")
rays = rc.ray_bundle_2d(
    [[0, 0, 0], [0, 1, 0], [0, 2, 0], [0, -1, 0], [0, -2, 0]], [0, 0.006, 1]
)
oec.ray_trace_2d(rays, [s102], out)

# Finding the paraxial focus by using a ray 0.1 mm above the optical
# axis and parallel to it.
r107 = rc.Ray([0, 0.1, 0], [0, 0, 1])
s103 = oec.SphericalRefraction(100, 0.03, 1, 1.5, 0.125)
oec.ray_trace_2d([r107], [s103], out101, focal_point=True)

# %% Task 11
""" 
Task 11: Testing the case of a point light source; theoretically all light 
rays (regardless of direction) should converge to the same point after 
propogating through the spherical lens.
"""
# Testing three rays originating from the same point but with different
# directions
r111 = rc.Ray([0, 2, 50], [0, 0, 1])
r112 = rc.Ray([0, 2, 50], [0, 0.01, 1])
r113 = rc.Ray([0, 2, 50], [0, -0.01, 1])
r114 = rc.Ray([0, 2, 50], [0, -0.02, 1])
r115 = rc.Ray([0, 2, 50], [0, -0.03, 1])
r116 = rc.Ray([0, 2, 50], [0, -0.04, 1])
s111 = oec.SphericalRefraction(100, 0.1, 1, 1.5, 3)


# Output plane at the position of paraxial focus (z=200mm).
out111 = oec.OutputPlane(155)

oec.ray_trace_2d([r111, r112, r113, r114, r115, r116], [s111], out111)
# %% Task 12
"""
Task 12: Trace a bundle of rays for a uniform collimated beam, 5mm to the
paraxial focal plane.
"""
out121 = oec.OutputPlane(200)
s121 = oec.SphericalRefraction(100, 0.03, 1, 1.5, 7.5)
s122 = oec.SphericalRefraction(95, 0, 1, 1.3, 10)

# Plotting the path of a collimated beam with 5 rings
# and a radius of 5mm through a refracting lens.

rays121 = rc.collimated_beam([0, 0, 0], 5, 5, [0, 0, 1])

oec.ray_trace_2d(rays121, [s121], out121)

# %% Task 13 and task 14
"""
Task 13 and task 14: For the same initial parameters above, plotting the spot 
diagram ofthe bundle of rays at the output plane, by using the optional 'spot' 
argument in the ray_trace_2d function. 

Used the optional 'rms' argument to find the RMS radius of the points on the
output plane from the optical axis.
"""
out131 = oec.OutputPlane(200)
s131 = oec.SphericalRefraction(100, 0.03, 1, 1.5, 7)

rays131 = rc.collimated_beam([0, 0, 0], 5, 5, [0, 0, 1])
oec.ray_trace_2d(rays131, [s131], out131, spot=True, rms=True)
plt.show()

# Cross-sectional view of the ray beam to check that it is uniform.
out132 = oec.OutputPlane(50)
rays132 = rc.collimated_beam([0, 0, 0], 5, 10, [0, 0, 1])
oec.ray_trace_2d(rays132, [], out131, spot=True)

# %% Task 15
"""
Part 1: Testing the 'arc_z' method which returns the z position of the maximum
and minimum point of the arc of a refracting surface.
"""

r151 = rc.Ray([0, 0, 0], [0, 0, 1])
s152 = oec.SphericalRefraction(100, 0.02, 1, 1, 10)
print(
    "Max aperture z position for this curved surface and is:%s mm" % (s152.arc_z(r151))
)

# %% Task 15
"""
Part 2: Testing the 'angle_find' method which uses the previous method to find
the angle subtended by the arc of the spherical aperture.

Also testing the 'plot' method which should plot the 
"""

print(
    "Angle substended by arc is from %s rad to %s rad"
    % (s152.angle_find(r151)[0], s152.angle_find(r151)[1])
)

# Plotting the arc.
s152.plot(r151)

# %% Task 15
"""
Part 3: Testing a plano-convex lens to see with a curvature of 0.02 mm^-1 and
a plane glass surface where the seperation of surfaces on the optical axis is
5 mm. Refractive index of the glass is 1.5168.

First testing with the plane surface facing the input.
"""
# Flat refracting surface at z = 95 mm.
s153 = oec.SphericalRefraction(95, 0, 1, 1.5168, 0.2)
# Spherical refracting surface with intercept at z = 100 mm.
s154 = oec.SphericalRefraction(100, -0.02, 1.5168, 1, 0.2)


# Finding an estimate for the paraxial focus using a ray 0.1 mm above the
# optical axis.
r152 = rc.Ray([0, 0.1, 80], [0, 0, 1])
out151 = oec.OutputPlane(250)
oec.ray_trace_2d([r152], [s153, s154], out151, focal_point=True)

# Now using a collimated beam with a radius of 10 mm and increasing the
# aperture radius.
rays151 = rc.collimated_beam([0, 0, 75], 10, 10, [0, 0, 1])
s153.__init__(95, 0, 1, 1.5168, 15)
s154.__init__(100, -0.02, 1.5168, 1, 15)
out151 = oec.OutputPlane(196.7487808269458)
oec.ray_trace_2d(rays151, [s153, s154], out151, spot=True, rms=True)
plt.show()

# Testing with a 5 mm collimated beam.
rays151 = rc.collimated_beam([0, 0, 75], 5, 10, [0, 0, 1])
s153.__init__(95, 0, 1, 1.5168, 7.5)
s154.__init__(100, -0.02, 1.5168, 1, 7.5)
oec.ray_trace_2d(rays151, [s153, s154], out151, spot=True, rms=True)
plt.show()

# Testing with a 2 mm collimated beam.
rays151 = rc.collimated_beam([0, 0, 75], 2, 10, [0, 0, 1])
s153.__init__(95, 0, 1, 1.5168, 3)
s154.__init__(100, -0.02, 1.5168, 1, 3)
oec.ray_trace_2d(rays151, [s153, s154], out151, spot=True, rms=True)
plt.show()

# Testing with a 0.1 mm collimated beam.
rays151 = rc.collimated_beam([0, 0, 75], 0.1, 10, [0, 0, 1])
s153.__init__(95, 0, 1, 1.5168, 0.15)
s154.__init__(100, -0.02, 1.5168, 1, 0.15)
oec.ray_trace_2d(rays151, [s153, s154], out151, spot=True, rms=True)
plt.show()

# %% Task 15
"""
Part 4: Repeating the same test above but now with the convex lens facing the
incident ray and the plane surface positioned after the convex lens.
"""
# Spherical refracting surface with intercept at z = 95 mm.
s155 = oec.SphericalRefraction(100, 0.02, 1, 1.5168, 0.2)
# Flat refracting surface with intercept at z = 100 mm.
s156 = oec.SphericalRefraction(105, 0, 1.5168, 1, 0.2)

# Paraxial focus estimate with 0.1 mm ray.
r153 = rc.Ray([0, 0.1, 80], [0, 0, 1])
out152 = oec.OutputPlane(250)
oec.ray_trace_2d([r153], [s155, s156], out152, focal_point=True)
plt.show()

# Testing with a 20 mm collimated beam.
rays152 = rc.collimated_beam([0, 0, 75], 20, 10, [0, 0, 1])
s155.__init__(100, 0.02, 1, 1.5168, 30)
s156.__init__(105, 0, 1.5168, 1, 30)
out152 = oec.OutputPlane(198.45270017751582)
oec.ray_trace_2d(rays152, [s155, s156], out152, spot=True, rms=True)
plt.show()

# Testing with a 10 mm collimated beam.
rays152 = rc.collimated_beam([0, 0, 75], 10, 10, [0, 0, 1])
s155.__init__(100, 0.02, 1, 1.5168, 15)
s156.__init__(105, 0, 1.5168, 1, 15)
out152 = oec.OutputPlane(198.45270017751582)
oec.ray_trace_2d(rays152, [s155, s156], out152, spot=True, rms=True)
plt.show()

# Testing with a 5 mm collimated beam.
rays152 = rc.collimated_beam([0, 0, 75], 5, 10, [0, 0, 1])
s155.__init__(100, 0.02, 1, 1.5168, 7.5)
s156.__init__(105, 0, 1.5168, 1, 7.5)
oec.ray_trace_2d(rays152, [s155, s156], out152, spot=True, rms=True)
plt.show()

# Testing with a 2 mm collimated beam.
rays152 = rc.collimated_beam([0, 0, 75], 2, 10, [0, 0, 1])
s155.__init__(100, 0.02, 1, 1.5168, 3)
s156.__init__(105, 0, 1.5168, 1, 3)
oec.ray_trace_2d(rays152, [s155, s156], out152, spot=True, rms=True)
plt.show()

# Testing with a 0.1 mm collimated beam.
rays152 = rc.collimated_beam([0, 0, 75], 0.1, 10, [0, 0, 1])
s155.__init__(100, 0.02, 1, 1.5168, 0.15)
s156.__init__(105, 0, 1.5168, 1, 0.15)
oec.ray_trace_2d(rays152, [s155, s156], out152, spot=True, rms=True)
plt.show()

# Testing with a 0.01 mm collimated beam.
rays152 = rc.collimated_beam([0, 0, 75], 0.01, 10, [0, 0, 1])
s155.__init__(100, 0.02, 1, 1.5168, 0.015)
s156.__init__(105, 0, 1.5168, 1, 0.015)
oec.ray_trace_2d(rays152, [s155, s156], out152, spot=True, rms=True)
plt.show()

# %% Lens Optimization Section
"""
Testing the 'lens_optimisation' module that contains s function to optimise 
the curvatures of a two surface lens system to minimise RMS.
"""
# Using a 5 mm radius (10 mm diamteter) collimated beam with the same density
# as for the 5 mm collimated beam in the section above.
# Setting the focal point to that of the plano convex lens in the section
# above.
rays161 = rc.collimated_beam([0, 0, 75], 5, 10, [0, 0, 1])
optimised_RMS = lop.lens_optimise(
    rays161, 100, 5, 1.5168, 10, [198.45270017751582, 198.45270017751582]
)
print(optimised_RMS)

# %% 10 mm beam
rays162 = rc.collimated_beam([0, 0, 75], 10, 10, [0, 0, 1])
optimised_RMS = lop.lens_optimise(
    rays162, 100, 5, 1.5168, 15, [198.45270017751582, 198.45270017751582]
)
print(optimised_RMS)

# %% 2 mm beam
rays163 = rc.collimated_beam([0, 0, 75], 2, 10, [0, 0, 1])
optimised_RMS = lop.lens_optimise(
    rays163, 100, 5, 1.5168, 3, [198.45270017751582, 198.45270017751582]
)
print(optimised_RMS)

# %% 0.1 mm beam
rays164 = rc.collimated_beam([0, 0, 75], 0.1, 10, [0, 0, 1])
optimised_RMS = lop.lens_optimise(
    rays164, 100, 5, 1.5168, 0.15, [198.45270017751582, 198.45270017751582]
)
print(optimised_RMS)

# %% 0.01 mm beam
rays165 = rc.collimated_beam([0, 0, 75], 0.01, 10, [0, 0, 1])
optimised_RMS = lop.lens_optimise(
    rays165, 100, 5, 1.5168, 0.15, [198.45270017751582, 198.45270017751582]
)
print(optimised_RMS)

# %% 20 mm beam
rays166 = rc.collimated_beam([0, 0, 75], 20, 10, [0, 0, 1])
optimised_RMS = lop.lens_optimise(
    rays166, 100, 5, 1.5168, 30, [198.45270017751582, 198.45270017751582]
)
print(optimised_RMS)


# %% Graphs
"""
Graphing the RMS for the plano-convex lens and the optimised best-form RMS
at different ray beam radii.
"""
Ray_Radii = [20, 10, 5, 2, 0.1, 0.01]
# RMS given to 3 s.f.
RMS_Plano_Convex = [0.596, 0.0683, 0.00836, 0.000531, 3.19e-8, 8.53e-9]
RMS_Optimised = [0.551, 0.0637, 0.00782, 0.000497, 3.16e-8, 8.46e-9]

plt.plot(Ray_Radii, RMS_Plano_Convex, label="Plano-convex", color="blue")
plt.plot(Ray_Radii, RMS_Optimised, label="Optimised", color="red")
plt.xlabel("Ray radius / mm")
plt.ylabel("RMS radius / mm")
plt.legend()
plt.grid()
plt.show()

# %% Graphs
"""Graphing the percentage decrease of RMS for the beam diameters."""
Percent_Change = [8.33, 6.60, 6.46, 6.40, 0.95, 0.82]

plt.plot(Ray_Radii, Percent_Change, "x", color="blue")
plt.xlabel("Ray radius / mm")
plt.ylabel("Percentage change")
plt.legend()
plt.grid()
plt.show()

# %% Thorlabs lens comparison.
"""
Comparing to Thorlab's LBF254-075-C lens.

Specifications
--------------
Radius of Curvature 1: 44.5 mm 
Radius of Curvature 2: -289.0 mm
Centre Thickness: 5.0 mm
Focal length: 75.0 mm
"""

# Finding the true focal position.
r161 = rc.Ray([0, 0.0001, 0], [0, 0, 1])
s_thor_1 = oec.SphericalRefraction(100, 1 / 44.5, 1, 1.5167, 7.5)
s_thor_2 = oec.SphericalRefraction(105, -1 / 289, 1.5167, 1, 7.5)
out161 = oec.OutputPlane(177.14)
oec.ray_trace_2d([r161], [s_thor_1, s_thor_2], out161, focal_point=True)

# RMS spread and spot diagram for the thor labs lens
rays_thor = rc.collimated_beam([0, 0, 0], 5, 10, [0, 0, 1])
oec.ray_trace_2d(rays_thor, [s_thor_1, s_thor_2], out161, rms=True, spot=True)

# %% Thorlabs lens comparison
"""
Finding the optimised lens at the paraxial focus of 177.14 mm and comparing
it with the Thorlabs lens.
"""
rays167 = rc.collimated_beam([0, 0, 75], 5, 10, [0, 0, 1])
# Constrained to a focal length of 75 mm
optimised_RMS = lop.lens_optimise(rays167, 100, 5, 1.5167, 7.5, [177.14, 177.14])
print(optimised_RMS)
