# %%
import pyvista as pv
from pyvista import examples
import numpy as np

# %%
r = 1
sphere = pv.Sphere(radius=r, center=(0, 0, 0), direction=(0, 0, 1), 
                   theta_resolution=30, phi_resolution=30, start_theta=0, 
                   end_theta=360, start_phi=0, end_phi=180)
outer_shell = pv.Tetrahedron(radius=r, center=(0, 0, 0))
pt = outer_shell.points[3]
cyl = pv.Cylinder(center=(0, 0, 0), direction=(0, 0, 1), 
                  radius=0.5*r, height=2*r, resolution=6, capping=True)
origin = np.array([0, 0, 0])
line_1 = pv.Line(origin, pt)
pl = pv.Plotter()
pl.add_mesh(sphere, show_edges=None, color=None, opacity=0.3)
pl.add_points(outer_shell.points, color='red',
              point_size=10)
pl.add_points(origin, color='black',
              point_size=15)
pl.add_mesh(line_1, color='w')
pl.add_mesh(cyl, style='wireframe', color='blue', line_width=2)

# rotate cylinder normal to sphere surface
import math
# angle_x = 90 - math.degrees(math.atan(pt[2]/pt[1]))
# angle_y = 90 - math.degrees(math.atan(pt[2]/pt[0]))
# angle_z = 90 - math.degrees(math.atan(pt[1]/pt[0]))
theta = math.degrees(math.acos(pt[2]/1))

cyl_rot = cyl.rotate_vector((pt[1], -pt[0], 0), -theta, point=cyl.center)
pl.add_mesh(cyl_rot, style='wireframe', color='red', line_width=2)

# cyl_rot_1 = cyl.rotate_x(-angle_x, point=cyl.center)
# pl.add_mesh(cyl_rot_1, style='wireframe', color='green', line_width=2)
# cyl_rot_2 = cyl_rot_1.rotate_y(angle_y, point=cyl.center)
# pl.add_mesh(cyl_rot_2, style='wireframe', color='red', line_width = 3)
# cyl_rot_y = cyl_rot_x.rotate_y(angle_y, point=cyl.center)
# pl.add_mesh(cyl_rot_y, style='wireframe', color='red', line_width = 5)
pl.show()


# %%
# translate to first tetrahedron point
pt = outer_shell.points[0]
v_tran = pt - cyl.center
cyl_tran = cyl.translate(v_tran, inplace=False)
pl.add_mesh(cyl_tran, style='wireframe', color='red')
pl.show()
# %%
# rotate cylinder normal to sphere surface
import math
angle_x = math.degrees(math.atan(pt[2]/pt[1]))
angle_y = math.degrees(math.atan(pt[2]/pt[0]))
angle_z = math.degrees(math.atan(pt[1]/pt[0]))

#%%
cyl_rot_x = cyl_tran.rotate_x(-angle_x, point=cyl_tran.center)
pl.add_mesh(cyl_rot_x, style='wireframe', color='green', line_width=2)
pl.show()
# %%
cyl_rot_y = cyl_rot_x.rotate_y(angle_y, point=cyl_tran.center)
pl.add_mesh(cyl_rot_y, style='wireframe', color='black', line_width = 5)
pl.show()
# %%
cyl_rot_z = cyl_rot_y.rotate_z(angle_z, point=cyl_tran.center)
pl.add_mesh(cyl_rot_z, style='wireframe', color='orange', line_width = 3)
pl.show()
# %%
# Test random perturbation of bullets
# 3/10/2023

# create rosette with one bullet

# Import packages
import pyvista as pv
import numpy as np
from math import pi
import trame
import math

# set geoemetric parameters
a = 0.3 # half max length across basal face
c =  1.5 # half max length across prism face
r0 = 1.0 # radius of center sphere
h0 = 0.25 # penetration depth of bullets
hp = 0.7 # heights of pyramid of bullets
n_arms = 1 # number of bullet arms

# set render parameters 
# pv.global_theme.restore_defaults()
bg_color = 'black' # background color of render
obj_color = 'white' # color of object
op = 0.9 # opacity of object

# create sphere 
sphere = pv.Sphere(radius=r0, center=(0, 0, 0), direction=(0, 0, 1), 
                   theta_resolution=30, phi_resolution=30, start_theta=0, 
                   end_theta=360, start_phi=0, end_phi=180)
r_outer = hp/2 + c - h0 + r0
outer_sphere = pv.Sphere(radius=r_outer, center=(0, 0, 0), direction=(0, 0, 1), 
                   theta_resolution=30, phi_resolution=30, start_theta=0, 
                   end_theta=360, start_phi=0, end_phi=180)
print(type(outer_sphere))

# create outer shell to "place" bullets on
# Modified fibbonaci lattice 
# Source: http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
epsilon = 0.33
goldenRatio = (1 + 5**0.5)/2
i = np.arange(0, n_arms) 
theta = 2 *pi * i / goldenRatio
phi = np.arccos(1 - 2*(i+epsilon)/(n_arms-1+2*epsilon))
x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)
outer_coords = r_outer*(np.column_stack((x, y, z)))

# create bullet arm
cyl = pv.Cylinder(center=(0.0, 0.0, c+hp), direction=(0.0, 0.0, -1.0), 
                  radius=2*a, height=2*c, resolution=6, capping=True)
pyr = pv.Cone(center=(0.0, 0.0, hp/2), direction=(0.0, 0.0, -1.0), 
              height=hp, radius=2*a, capping=True, angle=None, resolution=6)
cyl = cyl.triangulate()
pyr = pyr.triangulate()
bullet = cyl.boolean_union(pyr)

# copy, translate, and rotate bullets
origin = np.array([0, 0, 0])
pl = pv.Plotter()
pl.background_color = bg_color
pl.add_mesh(sphere, show_edges=None, color = obj_color, opacity=op)
# for i in range(len(outer_shell.points)):
for i in range(len(outer_coords)):
    # pt = outer_shell.points[i]
    pt = outer_coords[i]
    print(pt)
    # translate
    translate_vector = pt - bullet.center
    bullet_translated = bullet.translate(translate_vector, inplace=False)
    # rotate 
    theta = math.degrees(math.acos(pt[2]/r_outer))
    if (pt[0]==0 and pt[1]==0):
        bullet_final = bullet_translated.rotate_vector((0, pt[2], -pt[1]), -theta, point=bullet_translated.center)
    else:
        bullet_final = bullet_translated.rotate_vector((pt[1], -pt[0], 0), -theta, point=bullet_translated.center)
    pl.add_mesh(bullet_final, show_edges=None, color = obj_color, opacity=op)
pl.show()

# %%
# use rotations about axes
theta_x = np.random.randint(0, 360)
theta_y = np.random.randint(0, 360)
theta_z =np.random.randint(0, 360)

# rotate about x-axis
bullet_rotate_x = bullet_final.rotate_x(theta_x, inplace=False)
pl.add_mesh(bullet_rotate_x, show_edges=None, color = 'blue', opacity=op)

# rotate about y-axis
bullet_rotate_y = bullet_rotate_x.rotate_y(theta_y, inplace=False)
pl.add_mesh(bullet_rotate_y, show_edges=None, color = 'red', opacity=op)

# rotate about z-axis
bullet_rotate_z = bullet_rotate_y.rotate_z(theta_z, inplace=False)
pl.add_mesh(bullet_rotate_z, show_edges=None, color = 'green', opacity=op)

# for i in range(6):
#     if i > 0:
#         rot = bullet_rotate_z.rotate_x(60*i, inplace=False)
#         pl.add_mesh(rot, color='green')

# pl.show()


# %%
