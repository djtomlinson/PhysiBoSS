#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 10:44:08 2024

@author: D.J. Tomlinson
"""

import numpy as np
import csv
import matplotlib.pyplot as plt
from math import pi, cos, sin, asin

def generate_points_in_circle(radius, distance, randomness):
    # Calculate the number of points in x and y direction
    nx = int(radius / distance)
    ny = int(radius / distance)
    print("ny = %d" % ny, "nx = %d" % nx)

    # Generate the points
    points = []
    for i in range(-nx, nx + 1):
        for j in range(-ny, ny + 1):
            x = i * distance + np.random.uniform(-randomness, randomness)
            y = j * distance + (i % 2) * distance / 2 + np.random.uniform(-randomness, randomness)
            if np.sqrt(x**2 + y**2) <= radius:
                points.append((x, y, 0))  # z = 0, case = 0

    return points

# "diffusers, as you can set these particles to case = 1 and have case = 1 secrete chemicals at the edge of colony,
# it is preferred to maintain signal concetrations within the c++ code, but this can be a quick implementation
def generate_external_diffusers(distance, xMin, xMax, yMin, yMax, points=[]):
    # Place points on boundary to act as diffusers
    width = xMax - xMin
    height = yMax - yMin
    Nx = int(width / distance)
    Ny = int(height / distance)
    
    for k in range(0, Nx):
        x = k * distance + xMin
        points.append((x, yMin, 1)) #case = 1
        points.append((x, yMax, 1))
        
    for l in range(0, Ny):
        y = l * distance + yMin
        points.append((xMin, y, 1)) #case = 1
        points.append((xMax, y, 1))
    
    return points

def generate_edge_diffusers(radius, distance, thickness=1.0,ePoints=[]):
    # Place points on boundary to act as diffusers
    dTheta = asin(distance/(radius+thickness*distance))
    nTheta = int(2*pi/dTheta)
    for m in range(0, nTheta+1):
        x = (radius+thickness*distance) * cos(m*dTheta)
        y = (radius+thickness*distance) * sin(m*dTheta)
        ePoints.append((x,y, 1)) #case = 1
    return ePoints

# Test the function
frac = 0.75 #size as fraction of box size
xmin, xmax = -500, 500
ymin, ymax = -500, 500
radius = frac * np.min([xmax-xmin,ymax-ymin]) * 0.5 #as fraction of box size
distance = (6.25/462)*radius*2.0 #particle size / spacing between cells
randomness = 0
points = generate_points_in_circle(radius, distance, randomness)
#pointsWdiffusers = generate_external_diffusers(distance, xmin, xmax, ymin, ymax, points=points)
pointsWbounds = generate_edge_diffusers(radius, distance, thickness=1.0, ePoints=points)

# Write the points to a CSV file
with open('cells.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(pointsWbounds)

print("Points have been written to 'cells.csv'")
print("no. points = %d" % int(float(np.size(points))/4))

xPlot = []
yPlot = []
zPlot = []

#plot to check the points are what's expected
for setP in pointsWbounds:
    xPlot.append(setP[0])
    yPlot.append(setP[1])
    zPlot.append(setP[2])
    
plt.scatter(xPlot,yPlot)
plt.show()