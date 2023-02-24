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
