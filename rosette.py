# %%
import pyvista as pv
import numpy as np
from math import pi
import trame
import math

class Rosette:
    """
    Class representing bullet rosette ice crystals
    """
    def __init__(self, a, c, r0, h0, hp, n_arms):
        # geoemetric parameters
        self.a = a # half max length across basal face
        self.c =  c # half max length across prism face
        self.r0 = r0 # radius of center sphere
        self.h0 = h0 # penetration depth of bullets
        self.hp = hp # heights of pyramid of bullets
        self.n_arms = n_arms # number of bullet arms

        # create sphere 
        sphere = pv.Sphere(radius=r0, center=(0, 0, 0), direction=(0, 0, 1), 
                        theta_resolution=30, phi_resolution=30, start_theta=0, 
                        end_theta=360, start_phi=0, end_phi=180)
        r_outer = hp/2 + c - h0 + r0

        # create outer shell to "place" bullets on
        if n_arms == 2: # line
            outer_shell = pv.Line(pointa=(-r_outer, 0.0, 0.0), 
                                pointb=(r_outer, 0.0, 0.0), resolution=1)
            outer_coords = outer_shell.points
        elif n_arms == 4: # tetrahedron
            outer_shell = pv.Tetrahedron(radius=r_outer, center=(0.0, 0.0, 0.0))
            outer_coords = outer_shell.points
        elif n_arms == 6: # octahedron
            outer_shell = pv.Octahedron(radius=r_outer, center=(0.0, 0.0, 0.0))
            outer_coords = outer_shell.points
        elif n_arms == 8: # cube
            # Note: this may not be the optimal solution for n=8, check later
            l  = (2*r_outer)/(3**(1/2))
            outer_shell = pv.Cube(center=(0.0, 0.0, 0.0), x_length=l, 
                                y_length=l, z_length=l)
            outer_coords = outer_shell.points
        else: 
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
        rosette = sphere
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
            # perform union (to create single mesh)
            rosette = rosette.boolean_union(bullet_final)
        
        self.model = rosette # final 3d mesh model 

    def copy(self):
        """
        create a new instance
        with the same data as this instance
        """
        return Rosette(self.a, self.c, self.r0, self.h0, self.hp, self.n_arms)

    def plot(self, bg_color='black', obj_color='white', op=0.9):
        """
        Interactive PyVista visualization
        """
        pl = pv.Plotter(off_screen=True, window_size=[720,720])
        pl.background_color = bg_color
        pl.add_mesh(self.model, show_edges=None, color = obj_color, opacity=op)
        return pl

    def random_rotate(self):
        """
        Rotate rosette in a random orientation
        """
        rotated = self.copy()
        deg_x = np.random.randint(1, 360)
        deg_y = np.random.randint(1, 360)
        deg_z = np.random.randint(1, 360)
        rotated_model = self.model.rotate_x(deg_x, inplace=False)
        rotated_model.rotate_y(deg_y, inplace=True)
        rotated_model.rotate_z(deg_z, inplace=True)
        rotated.model = rotated_model
        return rotated

    def render(self, cam): 
        """
        Render orthographic (parallel) projection
        """
        pass

# %%
# Test instantiation
# set geoemetric parameters
a = 0.3 # half max length across basal face
c =  1.5 # half max length across prism face
r0 = 1.0 # radius of center sphere
h0 = 0.25 # penetration depth of bullets
hp = 0.7 # heights of pyramid of bullets
n_arms = 6 # number of bullet arms

# set render parameters 
# pv.global_theme.restore_defaults()
bg_color = 'black' # background color of render
obj_color = 'white' # color of object
op = 0.9 # opacity of object

# plot 
test = Rosette(a, c, r0, h0, hp, n_arms)
pl = test.plot()
pl.show(screenshot='test_render.png', interactive=False, jupyter_backend='none')
# pl.save_graphic('test_render.svg')

# # # rotate original then plot
# test_rotated = test.random_rotate()
# pl_rotated = test_rotated.plot()
# pl_rotated.show(screenshot='test_rotated_render.png')

#%%
# Create test subsample of particles
r0 = 1
hp = 0.7 # constant for now
h0 = 0.3 # constant for now
a_list = np.arange(0.1, 0.5, 0.1)
c_list = np.arange(0.5, 2.5, 0.5)
n_arms_list = np.arange(1, 10, 1)

# folder to save images
save_path = '/Users/josephko/research/ice_renders/sample_rosettes'

for a in a_list:
    for c in c_list:
        for n_arms in n_arms_list:
            ros = Rosette(a, c, r0, h0, hp, n_arms)
            ros_rotated = ros.random_rotate()
            pl = ros_rotated.plot()
            file_name = f'ros_n{n_arms}_a{a}_c{c}.png'
            file_path = save_path + '/' + file_name
            pl.show(screenshot=file_path, interactive=False, jupyter_backend='none')

# %%
