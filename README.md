Equal Earth Projection
======================
The Equal Earth map projection is a ``matplotlib`` add-on:

>The Equal Earth map projection is a new equal-area pseudocylindrical
>projection for world maps. It is inspired by the widely used Robinson
>projection, but unlike the Robinson projection, retains the relative size
>of areas. The projection equations are simple to implement and fast to
>evaluate. Continental outlines are shown in a visually pleasing and
>balanced way.

Projection developed by Bojan Šavrič (@BojanSavric):
* https://doi.org/10.1080/13658816.2018.1504949
* https://www.researchgate.net/publication/326879978_The_Equal_Earth_map_projection


![Example](charts/Equal_Earth_Tissot.png)

![Equations](charts/equations.jpeg)

where:

* λ and φ are the longitude and the latitude, respectively
* θ is a parametric latitude
* A1 =  1.340264
* A2 = -0.081106
* A3 =  0.000893
* A4 =  0.003796

Usage
-----
Importing the module causes the Equal Earth projection to be registered with
Matplotlib so that it can be used when creating a subplot::

    >>>import matplotlib.pyplot as plt
    >>>import EqualEarth
    >>>longs = [-110, 100, 100, -110]
    >>>lats = [40, 40, -40, 40]
    >>>fig = plt.figure('Equal Earth Projection')
    >>>ax = fig.add_subplot(111, projection="equal_earth")
    >>>ax.plot(np.deg2rad(longs), np.deg2rad(lats))
    >>>plt.grid(True)
    >>>plt.show()

![Result](charts/result.png)

Note that all data must be in radians, so be sure to use ``np.deg2rad()``
before plotting with data in degrees.

Sources
-------
Based on code from:
* https://matplotlib.org/gallery/misc/custom_projection.html

as well as code from @mbostock:
* https://beta.observablehq.com/@mbostock/equal-earth-projection


Issues
------
* Does not accept data in degrees, so data must be converted to radians first.
* The figure facecolor gets overdrawn with the axes facecolor

@Author: Dan Neuman (@dan613)

