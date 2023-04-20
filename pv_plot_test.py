
import pyvista as pv
from pyvista import examples

mesh = examples.load_airplane()

plotter = pv.Plotter()    # instantiate the plotter
plotter.add_mesh(mesh)    # add a mesh to the scene
plotter.show()   
