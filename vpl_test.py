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
mesh = mesh.Mesh.from_file(stl_file)

vpl.mesh_plot(mesh)

file_path = active_dir + '/' + 'test.png'
vpl.show()