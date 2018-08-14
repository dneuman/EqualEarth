#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Coastlines using Equal Earth Projection
=======================================

Example map showing how to use the Equal Earth Projection, and drawing
Tissot Indicatrices of Deformation to show how this projection distorts
land masses.

Requirements
------------
shapefile (from pyshp) is required to read the map data. This is available
from Anaconda, but must be installed first::

    >>>conda install shapefile

"""

import shapefile
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
import numpy as np
import EqualEarth  # Automatically registers in matplotlib on import

def GetAxes(figname, **kwargs):
    """
    Return matplotlib axes with Equal Earth projection
    """
    plt.tight_layout()  # Maximizes use of space
    fig = plt.figure(figname)
    fig.clear()
    ax = fig.add_subplot(111, projection='equal_earth', **kwargs)
    plt.grid(True)
    return ax

def DrawShapes(ax, sf, **kwargs):
    if sf.shapeType == shapefile.POLYGON:
        for shape in sf.shapes():
            verts = np.deg2rad(shape.points)
            patch = patches.Polygon(verts, **kwargs)
            ax.add_patch(patch)
    elif sf.shapeType == shapefile.POLYLINE:
        for shape in sf.shapes():
            verts = np.deg2rad(shape.points)
            path = patches.mlines.Path(verts)
            patch = patches.PathPatch(path, **kwargs)
            ax.add_patch(patch)

def DrawEllipse(ax, ll, width_deg, resolution=50):
    long, lat = ll
    # Use a path instead of the regular Ellipse patch to improve resolution
    height = np.deg2rad(width_deg)/2.  # use as radius, not diameter
    width = height/np.cos(lat)
    t = np.linspace(0., 2.*np.pi, resolution)
    t = np.r_[t, [0]]  # append starting point to close path
    longs = width * np.cos(t) + long
    lats = height * np.sin(t) + lat
    verts = np.column_stack([longs, lats])
    patch = patches.Polygon(verts,
                        facecolor='r', alpha=.4, edgecolor='none', zorder=5.)
    ax.add_patch(patch)

def DrawTissot(ax, width=10., resolution=50):
    """
    Draw Tissot Indicatrices of Deformation over the map projection to show
    how the projection deforms equally-sized circles at various points
    on the map.

    Parameters
    ----------
    ax : axes (subplot) to draw on
    width : (float default 5.) width of circles in degrees of latitude
    resolution : (int default 50) Number of points in circle
    """
    degrees = 30
    for lat in range(-degrees, degrees+1, degrees):
        for long in range(-180, 181, degrees):
            DrawEllipse(ax, np.deg2rad([long, lat]), width, resolution)
    for lat in [-60, 60]:
        for long in range(-180, 181, 2*degrees):
            DrawEllipse(ax, np.deg2rad([long, lat]), width, resolution)
    for lat in [-90, 90]:
        DrawEllipse(ax, np.deg2rad([0, lat]), width, resolution)
    ax.set_title('Equal Earth Projection with\n'
                 'Tissot Indicatrices of Deformation')


yellow = '#FEFEE6'
blue = '#CEEAFD'
matplotlib.rcParams['figure.facecolor'] = 'w'
matplotlib.rcParams['axes.facecolor'] = blue
matplotlib.rcParams['axes.edgecolor'] = 'k'
matplotlib.rcParams['grid.color'] = 'k'
matplotlib.rcParams['grid.alpha'] = .15

paths = ['ne_110m_land/ne_110m_land',
         'ne_110m_coastline/ne_110m_coastline',
         'ne_110m_lakes/ne_110m_lakes']
#      land    coast   lakes
ecs = ['none', 'k',    'k']  # edgecolors
fcs = [yellow, 'none', blue] # facecolors
z = 0.

ax = GetAxes('Equal Earth Tissot')
for path, ec, fc in zip(paths, ecs, fcs):
    sf = shapefile.Reader(path)
    DrawShapes(ax, sf, lw=.25, edgecolor=ec, facecolor=fc, zorder=z)
    z += .1
DrawTissot(ax)

plt.show()
