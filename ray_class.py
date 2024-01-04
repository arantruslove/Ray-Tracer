# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 16:57:16 2022

@author: HP
"""

"""Module that contains the Ray class."""

# MODLUES
import numpy as np


# CLASSES
class Ray:
    """
    Class that defines a light ray and is initialised with position and
    direction parameters.
    """

    def __init__(self, p, k):
        """


        Parameters
        ----------
        p : array_like
            List that contains the coordinates of the starting points of the
            ray.
        k : array_like
            List that describes the initial direction of travel of the ray.

        Returns
        -------
        None.

        """
        self.position = [np.array(p)]  # Intitial position of the ray.
        self.direction = [np.array(k)]  # Initial direction of the ray.

    def p(self):
        """Returns the current (last) position of the ray."""
        return self.position[-1]

    def k(self):
        """Returns the current (last) normalised direction of the ray."""
        return (1 / np.linalg.norm(self.direction[-1])) * self.direction[-1]

    def append(self, p=[0, 0, 0], k=[0, 0, 0]):
        """Method that appends a new position and direction to the ray."""
        self.position.append(np.array(p))
        self.direction.append(np.array(k))

    def vertices(self):
        """Method that returns all the vertices of the ray."""
        vertices = []
        for i in range(len(self.position)):
            # Stores vertex coordinates in a list rather than a numpy array.
            vertices.append(self.position[i].tolist())
        return vertices

    def focal_z(self):
        """
        Method that ignores the movement of the ray in the x direction and
        models the ray as a line in the form of y = mz + c to find the z
        coordinate where it crosses the optical axis.

        This method is used to find an estimate of the paraxial focus point
        from a single ray travelling close to and parallel to the optical axis.
        """

        direction = self.k()
        position = self.p()
        # Gradient of the ray in the y-z plane.
        # Gradient = change in y / change in z .
        gradient = direction[1] / direction[2]

        # Finding the constant in y = mz + c by comparing with the equation:
        # y - y1 = m(z - z1) and rearranging.
        constant = position[1] - gradient * position[2]

        # z-coordinate of the intercept point.
        z = -constant / gradient
        return z


def ray_bundle_2d(starting_points, direction):
    """
    Function that creates a 2d bundle of rays with the same direction but
    different starting points.

    Paramaters
    ----------
    starting_points : array_like
        Nested list that contains the coordinates of the starting points of
        each ray.
    drection : array_like
        Array that describes the direction that the rays travel in.
    """
    rays = []
    for i in range(len(starting_points)):
        ray1 = Ray(starting_points[i], direction)
        rays.append(ray1)
    return rays


def collimated_beam(centre, radius, number, direction):
    """
    Function that creates a list of rays to form a uniform, circular
    collimated beam.

    Parameters
    ---------
    centre: array_like
        Coordinates of the centre of the collimated beam.
    radius: float
        Radius of the collimated beam.
    number: int
        The number of rings of points that the beam has.
    direction: array_like
        Array that describes the direction of travel of the rays in the
        collimated beam.
    """

    s = radius / (number)  # Spacing between each ring of points.
    rays = []

    # Loops over each ring of the collimated beam.
    for i in range(number + 1):
        # First ray is at the centre of the collimated ray.
        if i == 0:
            rays.append(Ray(centre, direction))

        # Relationship between angle between rays for each ring of rays.
        else:
            theta = np.pi / (3 * i)

            # Converts from plane polar to cartesian coordinates in order
            # to determine the position of the rays relative to the centre
            # of the collimated ray.
            for j in range(6 * i):
                p1 = [
                    centre[0] + i * s * np.cos(j * theta),
                    centre[1] + i * s * np.sin(j * theta),
                    centre[2],
                ]
                rays.append(Ray(p1, direction))
    return rays
