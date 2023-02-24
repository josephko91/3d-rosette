# %%
# Import packages
import pyvista as pv
import numpy as np
from math import pi
import trame
import math

# set geometric parameters
a = 0.25 # half max length across basal face
c =  1.0 # half max length across prism face
r0 = 1.0 # radius of center sphere
h0 = 0.75 # penetration depth of bullets
hp = 1.0 # heights of pyramid of bullets
n_arms = 8 # number of bullet arms

# create sphere 
sphere = pv.Sphere(radius=r0, center=(0, 0, 0), direction=(0, 0, 1), 
                   theta_resolution=30, phi_resolution=30, start_theta=0, 
                   end_theta=360, start_phi=0, end_phi=180)
r_temp = hp/2 + c - h0 + r0
print(r_temp)

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

# # plot
# pl = pv.Plotter()
# pl.add_mesh(sphere, show_edges=None, color='white', opacity=0.5)
# pl.add_mesh(outer_shell, style='wireframe', color='red')
# # pl.add_points(outer_shell.points, color='red',
#             #   point_size=5)
# # pl.show()

# create bullet arm
cyl = pv.Cylinder(center=(0.0, 0.0, c+hp), direction=(0.0, 0.0, 1.0), 
                  radius=2*a, height=2*c, resolution=6, capping=True)
pyr = pv.Cone(center=(0.0, 0.0, hp/2), direction=(0.0, 0.0, -1.0), 
              height=hp, radius=2*a, capping=True, angle=None, resolution=6)
cyl = cyl.triangulate()
pyr = pyr.triangulate()
bullet = cyl.boolean_union(pyr)
# plotter = pv.Plotter()
# _ = plotter.add_mesh(cyl, color='r', style='wireframe', line_width=1)
# _ = plotter.add_mesh(pyr, color='b', style='wireframe', line_width=1)
# _ = plotter.add_mesh(bullet, color = 'white')
# plotter.show()

# # %%
# # plot
# pl = pv.Plotter()
# pl.add_mesh(sphere, style='wireframe', show_edges=True, color='white')
# # pl.add_mesh(outer_shell, style='wireframe', color='red')
# pl.add_points(outer_shell.points, color='red',
#               point_size=10)
# pl.add_mesh(bullet, style='wireframe', color='blue')
# pl.show()

# #%%
# # TEST: translate to first tetrahedron point
# origin = np.array([0, 0, 0])
# pl = pv.Plotter()
# pl.add_mesh(sphere, show_edges=None, color='white', opacity=0.5)
# pl.add_points(outer_shell.points, color='red',
#                 point_size=10)
# i=4
# pt = outer_shell.points[i]
# # line_1 = pv.Line(origin, pt)
# translate_vector = pt - bullet.center
# bullet_translated = bullet.translate(translate_vector, inplace=False)
# # pl.add_mesh(bullet_translated, style='wireframe', color='red')
# # pl.add_mesh(line_1, color='w')
# # # pl.show()
# theta = math.degrees(math.acos(pt[2]/r_temp))
# print((pt[1], -pt[0], 0))
# print(pt)
# print(theta)
# # %%
# # rotate bullet
# theta = math.degrees(math.acos(pt[2]/r_temp))
# # bullet_final = bullet_translated.rotate_vector((pt[1], -pt[0], 0), -theta, point=bullet_translated.center)
# if (pt[0]==0 and pt[1]==0):
#     bullet_final = bullet_translated.rotate_vector((0, pt[2], -pt[1]), -theta, point=bullet_translated.center)
# else:
#     bullet_final = bullet_translated.rotate_vector((pt[1], -pt[0], 0), -theta, point=bullet_translated.center)

# pl.add_mesh(bullet_final, style='wireframe', color='green', line_width=2)
# pl.show()


# translate to first tetrahedron point
origin = np.array([0, 0, 0])
pl = pv.Plotter()
pl.add_mesh(sphere, show_edges=None, color='white', opacity=0.9)
# pl.add_points(outer_shell.points, color='red',
#                 point_size=10)
for i in range(len(outer_shell.points)):
    pt = outer_shell.points[i]
    # line_1 = pv.Line(origin, pt)
    translate_vector = pt - bullet.center
    bullet_translated = bullet.translate(translate_vector, inplace=False)
    # pl.add_mesh(bullet_translated, style='wireframe', color='red')
    # pl.add_mesh(line_1, color='w')
    # # pl.show()

    # rotate bullet
    theta = math.degrees(math.acos(pt[2]/r_temp))
    if (pt[0]==0 and pt[1]==0):
        bullet_final = bullet_translated.rotate_vector((0, pt[2], -pt[1]), -theta, point=bullet_translated.center)
    else:
        bullet_final = bullet_translated.rotate_vector((pt[1], -pt[0], 0), -theta, point=bullet_translated.center)
    # pl.add_mesh(bullet_final, style='wireframe', color='green', line_width=2)
    pl.add_mesh(bullet_final, show_edges=None, color='white', opacity=0.9)
pl.show()

# # %%
# # test: place one bullet into sphere
# bullet_1 = bullet.translate([0, 0, -(2*c+hp+r0-h0)], inplace=False)
# sphere = sphere.triangulate()
# test_union = sphere.boolean_union(bullet_1)
# plotter_2 = pv.Plotter()
# _ = plotter_2.add_mesh(bullet_1, color='r', style='wireframe', line_width=1)
# _ = plotter_2.add_mesh(sphere, color='b', style='wireframe', line_width=1)
# _ = plotter_2.add_mesh(test_union, color = 'white')
# plotter_2.show()
# # %%
# cyl.points

# # %%
# from pyvista import examples
# pl = pv.Plotter()
# pl.show_axes()
# mesh1 = examples.download_teapot()
# _ = pl.add_mesh(mesh1)
# mesh2 = mesh1.flip_normal([1.0, 1.0, 1.0], inplace=False)
# _ = pl.add_mesh(mesh2, style='wireframe')
# pl.show()
# # %%
# mesh = pv.Cube(center=(1, 0, 0))
# rot = mesh.rotate_x(30, inplace=False)
# pl = pv.Plotter()
# _ = pl.add_mesh(rot)
# _ = pl.add_mesh(mesh, style='wireframe', line_width=3)
# _ = pl.add_axes_at_origin()
# pl.show()
# # %%
# mesh = pv.Cube(center=(1, 0, 0))
# rot = mesh.rotate_x(30, point=(1,0,0), inplace=False)
# pl = pv.Plotter()
# _ = pl.add_mesh(rot)
# _ = pl.add_mesh(mesh, style='wireframe', line_width=3)
# _ = pl.add_axes_at_origin()
# pl.show()
# # %%

# %%
