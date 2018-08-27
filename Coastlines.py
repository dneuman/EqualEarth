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

@Author: Dan Neuman (@dan613)
@Version: 1.1
@Date: 27 Aug 2018

New in this version
-------------------
* The main routine changes map colors without changing global colors
* The DrawCoastlines routine has simplified color arguments. The water color
  comes from the axes color, and the coastline and lake edge colors are the
  same. The land color is specified with ``facecolor`` or ``fc``, and the
  coastline and lake edge color is specified with ``edgecolor`` or ``ec``.
"""

import shapefile
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import EqualEarth  # Automatically registers in matplotlib on import

def GetAxes(figname='Equal Earth Projection', show_labels=True,
            figprops=None, **kwargs):
    """
    Return matplotlib axes with Equal Earth projection

    Parameters
    ----------
    figname : int or str, optional, default: 'Equal Earth Projection'
        Name of figure. If a number is given, will default to
        'Figure {figname}'
    show_labels : bool, optional, default: True
        Option to turn off grid labels.
    figprops : dict, optional, default: None
        Dictionary of properties to send to the figure during creation.
    kwargs : optional
        Keyword arguments will be sent to the plt.Axes object during creation.
    """
    fig = plt.figure(figname, **figprops)
    fig.clear()
    ax = fig.add_subplot(111, projection='equal_earth', **kwargs)
    if not show_labels:
            ax.yaxis.set_ticklabels([])
            ax.xaxis.set_ticklabels([])
    plt.tight_layout()  # Maximizes use of space
    plt.grid(True)
    return ax

def DrawShapes(ax, sf, **kwargs):
    """
    Draw shapes from the supplied shapefile onto the supplied axes

    Parameters
    ----------
    ax : axes (subplot)
        axes to draw on
    sf : shapefile.Reader object
        The shapefile containing the shapes to draw
    kwargs : optional
        Keyword arguments to send to the patch object. This will generally
        be edge and face colors, line widths, alpha, etc.
    """
    if sf.shapeType == shapefile.POLYGON:
        for shape in sf.shapes():
            verts = shape.points
            patch = patches.Polygon(verts, **kwargs)
            ax.add_patch(patch)
    elif sf.shapeType == shapefile.POLYLINE:
        for shape in sf.shapes():
            verts = shape.points
            path = patches.mlines.Path(verts)
            patch = patches.PathPatch(path, **kwargs)
            ax.add_patch(patch)

def DrawEllipse(ax, ll, width_deg, resolution=50):
    """
    Draw an ellipse on the supplied axes. Technically, a circle is drawn (an
    ellipse with equal height and width), but this usually becomes an ellipse
    on the projection axes.

    Parameters
    ----------
    ax : axes (subplot)
        axes to draw on
    ll : tuple of floats
        longitude and latitude coordinates (in degrees) to draw the ellipse
    width_deg : float
        Width of ellipse in degrees
    resolution : int, optional, default: 50
        number of points to use in drawing the ellipse
    """
    long, lat = ll
    # Use a path instead of the regular Ellipse patch to improve resolution
    height = width_deg/2.  # use as radius, not diameter
    width = height/np.cos(np.deg2rad(lat))
    t = np.linspace(0., 2. * np.pi, resolution)
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
    ax : axes (subplot)
        axes to draw on
    width : float, optional, default: 5.
        width of circles in degrees of latitude
    resolution : int, optional, default: 50
        Number of points in circle
    """
    degrees = 30
    for lat in range(-degrees, degrees+1, degrees):
        for long in range(-180, 181, degrees):
            DrawEllipse(ax, [long, lat], width, resolution)
    for lat in [-60, 60]:
        for long in range(-180, 181, 2*degrees):
            DrawEllipse(ax, [long, lat], width, resolution)
    for lat in [-90, 90]:
        DrawEllipse(ax, [0, lat], width, resolution)
    ax.set_title('Equal Earth Projection with\n'
                 'Tissot Indicatrices of Deformation')

def DrawCoastlines(ax, paths=None, edgecolor='k', facecolor='#FEFEE6',
                   linewidth=.25, **kwargs):
    """
    Draw land masses, coastlines, and major lakes. Colors and linewidth
    can be supplied.

    Parameters
    ----------
    ax : axes
        axes to draw on
    paths : list of str, optional, default: None
        List of paths to map data, if they aren't in the default location. The
        paths may be fully-specified or relative, and must be in format:
            ['land path', 'coastline path', 'lake path']
    edgecolor, ec : color, optional, default: black
        Color for coastlines and lake edges. ``ec`` can be used as a shortcut.
    facecolor, fc : color, optional, default: yellow
        Color for land. ``fc`` can be used as a shortcut.
    linewidth, lw : float, optional, default: .25
        Line width of coastlines and lake edges.
    """

    # Set up colors, overriding defaults if shortcuts given
    bc = ax.get_facecolor()           # background color
    ec = kwargs.get('ec', edgecolor)  # edge color
    fc = kwargs.get('fc', facecolor)  # face color
    lw = kwargs.get('lw', linewidth)  # line width

    #        land   coast   lakes
    edges = ['none', ec,    ec]
    faces = [fc,    'none', bc]

    if not paths:
        paths = ['maps/ne_110m_land/ne_110m_land',
                 'maps/ne_110m_coastline/ne_110m_coastline',
                 'maps/ne_110m_lakes/ne_110m_lakes']
    z = 0.
    for path, f, e in zip(paths, faces, edges):
        sf = shapefile.Reader(path)
        DrawShapes(ax, sf, linewidth=lw, zorder=z,
                   edgecolor=e, facecolor=f)
        z += .1

if __name__ == '__main__':
    blue = '#CEEAFD'

    ax = GetAxes('Equal Earth Tissot', show_labels=False,
                 figprops={'figsize': (10., 6.),
                           'facecolor': 'w'},
                 facecolor=blue)
    ax.spines['geo'].set_color('k')
    ax.grid(color='k', alpha=.15)

    DrawCoastlines(ax)
    DrawTissot(ax)

    plt.show()
