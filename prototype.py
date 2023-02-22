# %%
# Import packages
import pyvista as pv
import numpy as np
from math import pi
import trame

# %%
# set geometric parameters
a = 0.25 # half max length across basal face
c =  1.0 # half max length across prism face
r0 = 1.0 # radius of center sphere
h0 = 0.25 # penetration depth of bullets
hp = 1.0 # heights of pyramid of bullets

# %%
# create sphere 
sphere = pv.Sphere(radius=r0, center=(0, 0, 0), direction=(0, 0, 1), 
                   theta_resolution=30, phi_resolution=30, start_theta=0, 
                   end_theta=360, start_phi=0, end_phi=180)
sphere.plot(smooth_shading=False, color='white')
# %%
# create bullet arm 
cyl = pv.Cylinder(center=(0.0, 0.0, 0.0), direction=(0.0, 0.0, 1.0), 
                  radius=2*a, height=2*c, resolution=6, capping=True)
pyr = pv.Cone(center=(0.0, 0.0, 2*c+hp/2), direction=(0.0, 0.0, 1.0), 
              height=hp, radius=2*a, capping=True, angle=None, resolution=6)
plotter = pv.Plotter()
_ = plotter.add_mesh(cyl, 'r')
_ = plotter.add_mesh(pyr, 'b')
plotter.show()
# %%
