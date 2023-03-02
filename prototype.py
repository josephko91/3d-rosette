# %%
# Import packages
import pyvista as pv
import numpy as np
from math import pi
import trame
import math

# set geoemetric parameters
a = 0.25 # half max length across basal face
c =  2.0 # half max length across prism face
r0 = 1.0 # radius of center sphere
h0 = 0.50 # penetration depth of bullets
hp = 0.75 # heights of pyramid of bullets
n_arms = 6 # number of bullet arms

# set render parameters 
# pv.global_theme.restore_defaults()
bg_color = 'black' # background color of render
obj_color = 'white' # color of object
op = 0.9 # opacity of object

# create sphere 
sphere = pv.Sphere(radius=r0, center=(0, 0, 0), direction=(0, 0, 1), 
                   theta_resolution=30, phi_resolution=30, start_theta=0, 
                   end_theta=360, start_phi=0, end_phi=180)
r_temp = hp/2 + c - h0 + r0
print(r_temp)

# create outer shell to "place" bullets on
if n_arms == 2: # line
    outer_shell = pv.Line(pointa=(-r_temp, 0.0, 0.0), 
                          pointb=(r_temp, 0.0, 0.0), resolution=1)
elif n_arms == 4: # tetrahedron
    outer_shell = pv.Tetrahedron(radius=r_temp, center=(0.0, 0.0, 0.0))
elif n_arms == 6: # octahedron
    outer_shell = pv.Octahedron(radius=r_temp, center=(0.0, 0.0, 0.0))
elif n_arms == 8: # cube
    l  = (2*r_temp)/(3**(1/2))
    outer_shell = pv.Cube(center=(0.0, 0.0, 0.0), x_length=l, 
                          y_length=l, z_length=l)
else: 
    pass

# create bullet arm
cyl = pv.Cylinder(center=(0.0, 0.0, c+hp), direction=(0.0, 0.0, -1.0), 
                  radius=2*a, height=2*c, resolution=6, capping=True)
pyr = pv.Cone(center=(0.0, 0.0, hp/2), direction=(0.0, 0.0, -1.0), 
              height=hp, radius=2*a, capping=True, angle=None, resolution=6)
cyl = cyl.triangulate()
pyr = pyr.triangulate()
bullet = cyl.boolean_union(pyr)

# pl = pv.Plotter()
# # pl.add_mesh(sphere, show_edges=None, color='white', opacity=0.9)
# # pl.add_mesh(cyl, style='wireframe', color='red')
# # pl.add_mesh(pyr, style='wireframe', color='blue')
# pl.add_mesh(bullet)
# pl.show()

# translate and rotate bullets
origin = np.array([0, 0, 0])
pl = pv.Plotter()
pl.background_color = bg_color
pl.add_mesh(sphere, show_edges=None, color = obj_color, opacity=op)
for i in range(len(outer_shell.points)):
    pt = outer_shell.points[i]
    # translate
    translate_vector = pt - bullet.center
    bullet_translated = bullet.translate(translate_vector, inplace=False)
    # rotate 
    theta = math.degrees(math.acos(pt[2]/r_temp))
    if (pt[0]==0 and pt[1]==0):
        bullet_final = bullet_translated.rotate_vector((0, pt[2], -pt[1]), -theta, point=bullet_translated.center)
    else:
        bullet_final = bullet_translated.rotate_vector((pt[1], -pt[0], 0), -theta, point=bullet_translated.center)
    pl.add_mesh(bullet_final, show_edges=None, color = obj_color, opacity=op)
pl.show()

# %%
