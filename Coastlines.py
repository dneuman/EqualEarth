#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 12 19:42:07 2018

@author: dan
"""

import shapefile
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib
import numpy as np
import EqualEarth

def GetAxes(figname, **kwargs):
    """
    Return matplotlib axes with Equal Earth projection
    """
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
        verts = shape.points
        verts.append(verts[0])
        rads = np.deg2rad(verts)
        codes = [Path.MOVETO] + \
                (len(verts)-2) * [Path.LINETO] + \
                [Path.CLOSEPOLY]
        path = Path(rads, codes)
        patch = patches.PathPatch(path, **kwargs)
        ax.add_patch(patch)


matplotlib.rcParams['axes.facecolor'] = 'lightcyan'
matplotlib.rcParams['axes.edgecolor'] = 'k'
matplotlib.rcParams['grid.color'] = 'k'

ax = GetAxes('Equal Earth')
sf = GetCoastlines()
DrawShapes(ax, sf, lw=.5, edgecolor='k', facecolor='lightyellow', zorder=0.)
plt.show()
