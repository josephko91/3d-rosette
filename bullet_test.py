# %%
# Import packages
import pyvista as pv
import numpy as np
from math import pi
import trame
import math

# set geometric parameters
a = 0.4 # half max length across basal face
c =  1.0 # half max length across prism face
hp = 0.75 # heights of pyramid of bullets

# create bullet arm
cyl = pv.Cylinder(center=(0.0, 0.0, c+hp), direction=(0.0, 0.0, -1.0), 
                  radius=2*a, height=2*c, resolution=6, capping=True)
pyr = pv.Cone(center=(0.0, 0.0, hp/2), direction=(0.0, 0.0, -1.0), 
              height=hp, radius=2*a, capping=True, angle=None, resolution=6)
cyl = cyl.triangulate()
pyr = pyr.triangulate()
bullet = cyl.boolean_union(pyr)
pl = pv.Plotter()
pl.add_mesh(cyl, color='r', style='wireframe', line_width=2)
pl.add_mesh(pyr, color='b', style='wireframe', line_width=2)
pl.add_mesh(bullet, color='white', show_edges=True, opacity=0.5)
pl.show()
# %%
