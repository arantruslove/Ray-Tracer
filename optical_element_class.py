# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 17:16:57 2022

@author: HP
"""
"""
Module that contains a class for the spherical refracting surface and an
output plane.
"""

# MODULES
import numpy as np
import matplotlib.pyplot as plt

# CLASSES


class OpticalElement:
    """Class that defines an optical element such as a spherical surface."""

    def propogate_ray(self, Ray):
        """Propogates a ray through the optical element."""

        if self.intercept(Ray) is not None:
            I = self.intercept(Ray)
            k = self.refracted_direction(Ray)
            # Appending the new position and direction to the ray.
            return Ray.append(I, k)

        else:
            return None


class SphericalRefraction(OpticalElement):
    """
    Class that defines a spherical refracting surface that inherits from
    from the OpticalElement class.
    """

    def __init__(self, z0, C, n1, n2, AR):
        """
        Initialises the spherical surface with the information required to
        describe its shape and refractive properties.


        Parameters
        ----------
        z0 : float
            The z-position at which the refracting surface intercepts the
            optical axis.
        C : float
            The curvature of the lens, equal to 1/radius.
        n1 : float
            Refractive index of the medium before the surface.
        n2 : float
            Refractive index of the medium after the surface.
        AR : float
            Aperture radius of the lens which is the radius in which incident
            light rays can pass through the aperture.
        """
        self.__zintercept = np.array([0, 0, z0])
        self.__curvature = C
        self.__n1 = n1  # Refractive index before the surface.
        self.__n2 = n2  # Refractive index after the surface.
        # Perpendicular distance from the z-axis to the top of the surface.
        self.__ApertureRadius = AR

    def centre(self):
        """
        Method that returns the coordinates of the centre of the spherical
        surface.
        """
        if self.__curvature == 0:
            raise Exception("No centre defined for a planar surface.")

        else:
            return self.__zintercept + np.array([0, 0, 1 / self.__curvature])

    def intercept(self, Ray):
        """
        Returns the coordinates of the intercept of the ray with the
        refracting surface.
        """
        k = Ray.k()  # Current direction of the ray.

        # For the case of a flat aperture (C = 0).
        if self.__curvature == 0:
            # Distance between ray starting point and plane.
            l = ((self.__zintercept.tolist())[2] - Ray.p()[2]) / k[2]
            I = Ray.p() + l * k  # Intercept of the ray with the plane.

        # For the case of a spherical aperture.
        else:
            centre = self.centre()
            # Vector representing the direction between the centre of the
            # sphere and the starting point of the ray.
            r = Ray.p() - centre
            R = 1 / self.__curvature  # Radius of the spherical surface.

            # From the quadratic formula: when there are no real solutions to
            # the equation (b**2-4ac < 0).
            if (np.dot(r, k)) ** 2 - ((np.linalg.norm(r)) ** 2 - R**2) < 0:
                return None

            # When the intercept coordinates are real.
            else:
                # Calculates the distance between the starting point of
                # the ray and the intercept point of the line with the
                # spherical surface (can give negative and positive
                # solutions).
                l1 = -np.dot(r, k) + np.sqrt(
                    (np.dot(r, k)) ** 2 - ((np.linalg.norm(r)) ** 2 - R**2)
                )
                l2 = -np.dot(r, k) - np.sqrt(
                    (np.dot(r, k)) ** 2 - ((np.linalg.norm(r)) ** 2 - R**2)
                )

                # Coordinates of the two theoretical intercept points.
                I1 = Ray.p() + l1 * k
                I2 = Ray.p() + l2 * k

                if l1 < 0:
                    # Two negative solutions means that the is
                    # ray travelling away from the spherical surface.
                    return None

                elif l2 < 0:
                    # Will pick I1 since it corresponds to l1 which is the
                    # only positive distance.
                    I = I1

                else:
                    # Both solutions are positive so will pick either the
                    # solution that meets the aperture radius criteria or
                    # if both solutions meet the criteria then will pick
                    # the smaller solution as this is the first intercept.
                    if abs(I2[1]) < self.__ApertureRadius:
                        I = I2

                    elif abs(I1[1]) < self.__ApertureRadius:
                        I = I1

                    else:
                        return None

        if abs(I[1]) < self.__ApertureRadius:
            return I
        else:
            return None

    def refracted_direction(self, Ray):
        """
        Determines the direction of the refracted ray and appends this
        direction and intercept to the ray.

        3d Snell's law equation used to find the direction of the ray
        after propogating through the aperture. Derivation of the 3d
        Snell's law equation found here:
        https://physics.stackexchange.com/questions/435512/snells-law-in-vector-form
        """
        if self.intercept(Ray) is not None:
            I = self.intercept(Ray)
            k1 = Ray.k()  # Normalised direction vector of incident ray.
            mu = self.__n1 / self.__n2  # Ratio of refractive indices.

            if self.__curvature == 0:
                # Vector normal to the x-y plane.
                n = np.array([0, 0, 1])
                test = np.dot(k1, n)

                # Case that normal to the plane and ray direction make an
                # acute angle.
                if test > 0:
                    pass

                # Need to reverse the direction of the normal when the
                # angle between the ray and the normal is an obtuse angle.
                else:
                    n = -n

            else:
                # Direction vector between center of sphere and intercept
                # point.
                Radial = self.centre() - I
                # Normalised vector that represents the normal to the
                # sphere at the intercept.
                n = -(1 / np.linalg.norm(Radial)) * Radial
                test = np.dot(k1, n)

                if test > 0:
                    pass
                else:
                    n = -n

            # Finding the angle between the direction of the incident ray
            # and the normal to the surface.
            theta = np.arccos(np.dot(n, k1))

            # Condition for TIR; when angle of incidence is greater than
            # the critical angle.
            if np.sin(theta) > 1 / mu:
                raise Exception("No propogation due to TIR.")
            else:
                # Unormalised ray direction vector after the surface.
                k2_unnormalised = np.sqrt(
                    1 - mu**2 * (1 - (np.dot(n, k1)) ** 2)
                ) * n + mu * (k1 - (np.dot(n, k1)) * n)
                # Normalised direction vector
                k2 = (1 / np.linalg.norm(k2_unnormalised)) * k2_unnormalised

                return k2
        else:
            return Ray.k()

    # Next three functions are solely for the purpose of plotting the graph
    # of a refracting surface and are not usually meant to be used alone.

    def arc_z(self, Ray):
        """
        Methods that returns the maximum and minimum height coordinates of the
        aperture.
        """
        r = self.__ApertureRadius
        R = abs(1 / self.__curvature)

        if r < R:
            # For the case of a curved surface.
            if self.__curvature != 0:
                # Possible z coordinates that satisfy the aperture radius.
                z = [
                    self.centre()[2] + np.sqrt(R**2 - r**2),
                    self.centre()[2] - np.sqrt(R**2 - r**2),
                ]

                # There will be two possible arcs of the sphere satisfying the
                # aperture radius, so a is used to find the arc that is
                # associated  with the intercept of the aperture with the
                # surface.
                a = (
                    abs(self.intercept(Ray)[2] - z[0]),
                    abs(self.intercept(Ray)[2] - z[1]),
                )

                if a[0] < a[1]:
                    zcoord = z[0]
                else:
                    zcoord = z[1]
                return zcoord

            # For the case of a flat surface.
            else:
                zcoord = self.__zintercept[2]
                return zcoord

        else:
            return None

    # Following three methods are solely meant for plotting the graph of the
    # refracting surfaces and are not usually meant to be accessed
    # seperately.

    def angle_find(self, Ray):
        """Finds the angles between which the arc of the aperture subtends."""
        zcoord = self.arc_z(Ray)
        r = self.__ApertureRadius
        centre = self.centre()

        # Distance between centre and z coordinate for geometric purposes.
        zdistance = abs(zcoord - centre[2])
        theta = np.arctan(r / (zdistance))

        # For the case where the aperture is ahead of the centre.
        if zcoord > centre[2]:
            return [-theta, theta]

        # When the aperture is behind the centre.
        else:
            return [np.pi - theta, np.pi + theta]

    def arc_plot(self, Ray):
        """
        Method that plots an arc for the case where curvature is not zero.
        """
        R = abs(1 / self.__curvature)
        r = self.__ApertureRadius
        centre = self.centre()
        test = self.intercept(Ray)[2] - centre[2]

        # When the aperture radius is greater than the radius of the
        # sphere.
        if r >= R and test > 0:
            # Angles representing a semi-circle when intercept is to the
            # right of the centre.
            angles = np.linspace(-0.5 * np.pi, 0.5 * np.pi)

        elif r >= R and test < 0:
            # When the intercept is behind the centre.
            angles = np.linspace(0.5 * np.pi, 1.5 * np.pi)

        else:
            # Finds an array of evenly spaced angles of the arc.
            thetas = self.angle_find(Ray)
            angles = np.linspace(thetas[0], thetas[1], 1000)

        # Uses the list of angles to find the coordinates located on the
        # arc.
        x = centre[2] + R * np.cos(angles)
        y = centre[1] + R * np.sin(angles)

        plt.plot(x, y, color="red", zorder=2)

    def plot(self, Ray):
        """Method that plots the aperture for all possible curvature."""

        # No plot when there is no intercept between ray and surface.
        if self.intercept(Ray) is None:
            return None

        else:
            z_intercept = self.__zintercept[2]
            r = self.__ApertureRadius
            # Straight, vertical line displayed for zero curvature.
            if self.__curvature == 0:
                plt.plot([z_intercept, z_intercept], [-r, r], color="red", zorder=2)

            else:
                self.arc_plot(Ray)


class OutputPlane(OpticalElement):
    """
    Output plane used to calculate a final point (vertex) for a ray which
    is useful for graphing the journey of the ray.
    """

    def __init__(self, z0):
        """

        Parameters
        ----------
        z0 : float
            z-coordinate of the position where the surface intercepts the
            optical axis.
        """
        self.zintercept = z0

    def intercept(self, Ray):
        """
        Returns the coordinates of the intercept of the ray with the output
        plane.
        """
        k = Ray.k()
        # Distance between ray starting point and plane.
        l = (self.zintercept - Ray.p()[2]) / k[2]
        I = Ray.p() + l * k  # Intercept of ray with plane.

        return I

    def end_ray(self, Ray):
        """
        Appends the coordinates of the intercept of the ray with the output
        plane to the ray.
        """
        I = self.intercept(Ray)
        k = Ray.k()
        return Ray.append(I, k)


# FUNCTIONS


def ray_trace_2d(
    rays,
    refracting_surfaces,
    output_plane,
    vertices=False,
    spot=False,
    rms=False,
    focal_point=False,
):
    """

    Function that inputs rays and surfaces and plots the path that the ray
    will take.


    Parameters
    ----------
    rays : array_like
        A list that contains rays.
    refracting_surfaces : array_like
        A list that contains refracting surfaces such as spherical surfaces.
    output_plane : object
        An object that is the plane in which the rays collide with at the end
        of their journey.
    vertices : bool, optional
        If true, lists all the vertices of each ray.
    spot : bool, optional
        If true, plots all the points of each ray on the output plane.
    rms: bool, optional
        If true, prints the value of the RMS distance of all the points from
        the optical axis at the output plane. Meant for rays that travel
        parallel to the optical axis.
    optical_intercept : bool,optionL
        Finds the intercept of a ray with the optical axis (when a single ray
        is input). Can be used to find the estimate of the position of the
        paraxial focus.
    """
    x_end = []
    y_end = []
    for i in range(len(rays)):
        # Propogates the rays through the apertures and stops the ray at the
        # output plane.
        for j in range(len(refracting_surfaces)):
            refracting_surfaces[j].plot(rays[i])
            refracting_surfaces[j].propogate_ray(rays[i])
        output_plane.end_ray(rays[i])

        # Extracting the y and z coordinates of the ray and placing them on
        # the y and x coordinates respectively on a graph.
        zcoordinates = []
        for j in range(len(rays[i].vertices())):
            zcoordinates.append(rays[i].vertices()[j][2])

        ycoordinates = []
        for j in range(len(rays[i].vertices())):
            ycoordinates.append(rays[i].vertices()[j][1])

            # Forms a list of coordinates for the x and y positions of the rays
            # at the output plane.
            x_end.append(rays[i].vertices()[-1][0])
            y_end.append(rays[i].vertices()[-1][1])

        # Plotting vertex coordinates on a graph.
        plt.ylabel("y / mm")
        plt.xlabel("z / mm")
        plt.errorbar(zcoordinates, ycoordinates, color="blue", zorder=1)
    plt.grid()
    plt.show()

    # CODE FOR THE OPTIONAL KEY WORD ARGUMENTS.

    if vertices is True:
        print("Vertices of ray %s:%s" % (i + 1, rays[i].vertices()))
    else:
        pass

    # Plotting the spot diagram at the output plane.
    if spot is True:
        plt.plot(x_end, y_end, ".")
        plt.xlabel("x-axis")
        plt.ylabel("y-axis")

        # Equal aspect ratio for the axes so that circular spot diagrams do
        # not appear elliptical.
        plt.gca().set_aspect("equal", adjustable="box")
        plt.draw()
    else:
        pass

    # Finding the RMS spot radius
    if rms is True:
        squares = 0
        for i in range(len(x_end)):
            # Sum of the distances between the points (x,y) and the optical
            # axis.
            squares += x_end[i] ** 2 + y_end[i] ** 2

        mean_squares = squares / len(x_end)  # Mean of the squares
        rms = np.sqrt(mean_squares)
        print("RMS radius from paraxial axis is: %s mm" % (rms))
        return rms
    else:
        pass

    # Finding the z-position of the paraxial focus.
    if focal_point is True and len(rays) == 1:
        z = rays[0].focal_z()
        print("z-position of the paraxial focus is %s mm" % (z))

    else:
        pass
