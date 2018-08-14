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
from matplotlib.path import Path
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

def GetCoastlines(path=None):
    """
    Return the shapefile containing coastlines
    """
    if not path:
        path = 'ne_110m_land/ne_110m_land'
    sf = shapefile.Reader(path)
    return sf

def DrawShapes(ax, sf, **kwargs):
    """
    Draw the shapes in the supplied shapefile
    """
    for shape in sf.shapes():
        if shape.shapeType != 5: continue
        verts = shape.points
        codes = [Path.MOVETO] + \
                (len(verts)-1) * [Path.LINETO]
        path = Path(np.deg2rad(verts), codes)
        patch = patches.PathPatch(path, **kwargs)
        ax.add_patch(patch)

def DrawEllipse(ax, ll, width):
    R = 6371.  # radius of Earth in km
    long, lat = ll
    # circumference is 2πR, so angle of longitude is:
    y = width/(2.*np.pi*R)
    x = y/np.cos(lat)
    p = patches.Ellipse((long, lat), x, y,
                        color='r', alpha=.4, edgecolor=None, zorder=5.)
    ax.add_patch(p)

def DrawTissot(ax):
    degrees = 30
    width = 7500.
    for lat in range(-degrees, degrees+1, degrees):
        for long in range(-180, 181, degrees):
            DrawEllipse(ax, np.deg2rad([long, lat]), width)
    for lat in [-60, 60]:
        for long in range(-180, 181, 2*degrees):
            DrawEllipse(ax, np.deg2rad([long, lat]), width)


matplotlib.rcParams['figure.facecolor'] = 'w'
matplotlib.rcParams['axes.facecolor'] = '#CEEAFD'
matplotlib.rcParams['axes.edgecolor'] = 'k'
matplotlib.rcParams['grid.color'] = 'k'
matplotlib.rcParams['grid.alpha'] = .15

ax = GetAxes('Equal Earth')
sf = GetCoastlines()
DrawShapes(ax, sf, lw=.5, edgecolor='k', facecolor='#FEFEE6', zorder=0.)
DrawTissot(ax)
ax.set_title('Equal Earth Projection with Tissot Indicatrices of Deformation')
plt.tight_layout()
plt.show()
