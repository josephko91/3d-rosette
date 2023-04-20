import pyvista as pv
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
import vtkplotlib as vpl
# Generate random views of the particle and save as png
active_dir = '/Users/josephko/research/ice_renders/20230413'
n_renders = 1_000

# Create a new plot
figure = pyplot.figure()
axes = figure.add_subplot(projection='3d')

# Load the STL files and add the vectors to the plot
stl_file = active_dir + '/' + 'ros_n8_a0.40_c1.00_9.stl'
your_mesh = mesh.Mesh.from_file(stl_file)
axes.add_collection3d(mplot3d.art3d.Poly3DCollection(your_mesh.vectors, color='black'))

# Auto scale to the mesh size
scale = your_mesh.points.flatten()
axes.auto_scale_xyz(scale, scale, scale)

pyplot.axis('off')

# Show the plot to the screen
file_path = active_dir + '/' + 'test.png'
pyplot.savefig(file_path)

# cube = pv.Cube()
# pl = pv.Plotter(off_screen=True)
# pl.add_mesh(cube)
# file_name = f'test.png'
# file_path = save_path + '/' + file_name
# image = pl.screenshot(None, return_img=True)
# print(type(image))
# # pl.show(auto_close=False)
# # image = pl.screenshot(None, return_img=True)
# # pl.close()

# print(type(image))

# for i in range(n_renders):
#     sphere = pv.Cube()
#     pl = pv.Plotter(off_screen=True)
#     pl.add_mesh(sphere)
#     file_name = f'test_{i}.png'
#     file_path = save_path + '/' + file_name
#     pl.screenshot(file_path, return_img=False)