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
    ll = np.deg2rad(ll)
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

def DrawCoastlines(ax, paths=None, edgedict=None, facedict=None,
                   linewidth=.25):
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
    edgedict : dict, optional, default: None
        Optional dictionary for line colors.
        The dict keys must be any of ['land', 'coast', 'lakes']
    facedict : dict, optional, default: None
        Optional dictionary for interior colors.
        The dict keys must be any of ['land', 'coast', 'lakes']
    linewidth : float, optional, default: .25
        Line width of coastlines
    """
    # Set up default colors
    yellow = '#FEFEE6' # wikipedia map colors (roughly)
    blue = '#CEEAFD'
    keys =               ['land', 'coast', 'lakes']
    ecd = dict(zip(keys, ['none', 'k',     'k']))  # edgecolors
    fcd = dict(zip(keys, [yellow, 'none',  blue])) # facecolors
    # Update with any supplied colors
    if edgedict: ecd.update(edgedict)
    if facedict: fcd.update(facedict)
    if not paths:
        paths = ['ne_110m_land/ne_110m_land',
                 'ne_110m_coastline/ne_110m_coastline',
                 'ne_110m_lakes/ne_110m_lakes']
    z = 0.
    for path, key in zip(paths, keys):
        sf = shapefile.Reader(path)
        DrawShapes(ax, sf, lw=linewidth, zorder=z,
                   edgecolor=ecd[key], facecolor=fcd[key])
        z += .1

if __name__ == '__main__':
    yellow = '#FEFEE6'
    blue = '#CEEAFD'

    matplotlib.rcParams['figure.facecolor'] = 'w'
    matplotlib.rcParams['axes.facecolor'] = blue
    matplotlib.rcParams['axes.edgecolor'] = 'k'
    matplotlib.rcParams['grid.color'] = 'k'
    matplotlib.rcParams['grid.alpha'] = .15

    ax = GetAxes('Equal Earth Tissot', show_labels=False,
                 figprops={'figsize': (10., 6.)})
    DrawCoastlines(ax)
    DrawTissot(ax)

    plt.show()
