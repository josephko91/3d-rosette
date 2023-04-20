import pyvista as pv
import numpy as np

# # Generate random views of a cube and save render as png
# save_path = '/Users/josephko/research/ice_renders/test'
# n = 100
# cube = pv.Cube()

# for i in range(n):
#     pl = pv.Plotter(off_screen=True)
#     deg_x = np.random.randint(1, 360)
#     deg_y = np.random.randint(1, 360)
#     deg_z = np.random.randint(1, 360)
#     rotated_cube = cube.rotate_x(deg_x, inplace=False)
#     rotated_cube.rotate_y(deg_y, inplace=True)
#     rotated_cube.rotate_z(deg_z, inplace=True)
#     pl.add_mesh(rotated_cube)
#     file_name = f'cube_{i}.png'
#     file_path = save_path + '/' + file_name
#     pl.screenshot(file_path, return_img=False)
#     pl.close()

# Generate random views of a cube and save render as png
save_path = '/Users/josephko/research/ice_renders/test'
n = 1_000
cube = pv.Cube()

pl = pv.Plotter(off_screen=True)
# Use name to replace mesh
pl.add_mesh(cube)

for i in range(n):
    deg_x = np.random.randint(1, 360)
    deg_y = np.random.randint(1, 360)
    deg_z = np.random.randint(1, 360)
    cube.rotate_x(deg_x, inplace=True)
    cube.rotate_y(deg_y, inplace=True)
    cube.rotate_z(deg_z, inplace=True)
    pl.add_mesh(cube)
    file_name = f'cube_{i}.png'
    file_path = save_path + '/' + file_name
    pl.screenshot(file_path, return_img=False)

pl.close()