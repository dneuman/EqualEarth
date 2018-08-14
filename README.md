Equal Earth Projection
======================

Based on code from:
* https://matplotlib.org/gallery/misc/custom_projection.html

and projection described by Bojan Šavrič (@BojanSavric):
* https://doi.org/10.1080/13658816.2018.1504949
* https://www.researchgate.net/publication/326879978_The_Equal_Earth_map_projection

as well as code from @mbostock:
* https://beta.observablehq.com/@mbostock/equal-earth-projection

Abstract:
    "The Equal Earth map projection is a new equal-area pseudocylindrical
     projection for world maps. It is inspired by the widely used Robinson
     projection, but unlike the Robinson projection, retains the relative size
     of areas. The projection equations are simple to implement and fast to
     evaluate. Continental outlines are shown in a visually pleasing and
     balanced way."

![Example](charts/Equal_Earth_Tissot.png)

Usage
-----
Importing the module causes the Equal Earth projection to be registered with
Matplotlib so that it can be used when creating a subplot::

    >>>import matplotlib.pyplot as plt
    >>>import EqualEarth
    >>>longs = [-110, 50, 50]
    >>>lats = [20, 20, -40]
    >>>fig = plt.figure('Equal Earth Projection')
    >>>ax = fig.add_subplot(111, projection="equal_earth")
    >>>ax.plot(np.deg2rad(longs), np.deg2rad(lats))
    >>>plt.show()

Note that all data must be in radians, so be sure to use ``np.deg2rad()``
before plotting with data in degrees.

Issues
------
* Does not accept data in degrees, so data must be converted to radians first.
* The figure facecolor gets overdrawn with the axes facecolor

@Author: Dan Neuman (@dan613)

