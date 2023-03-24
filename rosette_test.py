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
        self.c = c # half max length across prism face
        self.r0 = r0 # radius of center sphere
        self.h0 = h0 # penetration depth of bullets
        self.hp = hp # heights of pyramid of bullets
        self.n_arms = n_arms # number of bullet arms

        # create sphere 
        sphere = pv.Sphere(radius=r0, center=(0, 0, 0), direction=(0, 0, 1), 
                        theta_resolution=30, phi_resolution=30, start_theta=0, 
                        end_theta=360, start_phi=0, end_phi=180)
        self.sphere = sphere

        # create outer shell to "place" bullets on
        r_outer = hp/2 + c - h0 + r0
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
        self.bullets = {} # save bullets in nested dictionary
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
            
            # add bullet attributes and mesh to dictionary
            bullet_entry = {}
            bullet_entry['mesh'] = bullet_final
            bullet_entry['xy_scale_factor'] = 1.0
            bullet_entry['z_scale_factor'] = 1.0
            bullet_entry['anchor_point'] = pt
            self.bullets[i] = bullet_entry

    def unify_mesh(self): 
        """
        Create single mesh using bean union operation
        """
        rosette = self.sphere
        for i in range(self.n_arms):
            bullet = self.bullets[i]
            bullet_mesh = bullet['mesh']
            rosette = rosette.boolean_union(bullet_mesh)
        
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

        if hasattr(self, 'model'):
            pl.add_mesh(self.model, show_edges=None, color = obj_color, opacity=op)
        else: 
            pl.add_mesh(self.sphere, show_edges=None, color = obj_color, opacity=op)
            for i in range(self.n_arms):
                bullet = self.bullets[i]['mesh']
                pl.add_mesh(bullet, show_edges=None, color = obj_color, opacity=op)
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
    
    def randomize_bullets(self, scaling=False, location=False, inplace=True):
        """
        Randomly perturb the scaling and location of bullets
        """
        if inplace:
            rosette = self
        else:
            rosette = self.copy()

        # Perturb bullet scaling
        if scaling == False:
            pass
        else:
            sf_basal = np.random.uniform(scaling[0], scaling[1]) # basal scaling factor
            sf_prism = np.random.uniform(scaling[2], scaling[3]) # prism scaling factor
            r_outer = rosette.hp/2 + rosette.c - rosette.h0 + rosette.r0
            for i in range(rosette.n_arms):
                bullet = rosette.bullets[i]
                pt = bullet['anchor_point']
                theta = math.degrees(math.acos(pt[2]/r_outer))
                # rotate parallel to north unit vector
                if (pt[0]==0 and pt[1]==0):
                    bullet['mesh'].rotate_vector((0, pt[2], -pt[1]), theta, point=[0,0,0], inplace=True)
                else:
                    bullet['mesh'].rotate_vector((pt[1], -pt[0], 0), theta, point=[0,0,0], inplace=True)
                # scale 
                t_mat = np.array([[sf_basal, 0, 0, 0],
                                  [0,sf_basal, 0, 0],
                                  [0, 0, sf_prism, 0],
                                  [0, 0, 0, 1]])
                bullet['mesh'].transform(t_mat, inplace=True)
                # rotate back to place
                if (pt[0]==0 and pt[1]==0):
                    bullet['mesh'].rotate_vector((0, pt[2], -pt[1]), -theta, point=[0,0,0], inplace=True)
                else:
                    bullet['mesh'].rotate_vector((pt[1], -pt[0], 0), -theta, point=[0,0,0], inplace=True)

        # Perturb bullet location
        if location==False:
            pass
        else:
            pass

        return rosette

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
pl.show()
# pl.show(screenshot='test_render.png', interactive=False, jupyter_backend='none')
# pl.save_graphic('test_render.svg')

# # # rotate original then plot
# test_rotated = test.random_rotate()
# pl_rotated = test_rotated.plot()
# pl_rotated.show(screenshot='test_rotated_render.png')

#%%
# 3/23/23
# Testing for randomizing bullets 
# print(test.__dir__())
# for key in test.bullets:
#     bullet = test.bullets[key]['mesh']
#     bullet.plot(jupyter_backend='static')

test_bullet = test.bullets[5]
print(test_bullet['anchor_point'])
pt = test_bullet['anchor_point']
r_outer = hp/2 + c - h0 + r0
theta = math.degrees(math.acos(pt[2]/r_outer))
if (pt[0]==0 and pt[1]==0):
    test_bullet_rotated = test_bullet['mesh'].rotate_vector((0, pt[2], -pt[1]), theta, point=[0,0,0])
else:
    test_bullet_rotated = test_bullet['mesh'].rotate_vector((pt[1], -pt[0], 0), theta, point=[0,0,0])
pl2 = pv.Plotter()
pl2.add_mesh(test_bullet['mesh'], show_edges=None, color = obj_color, opacity=op)
pl2.add_mesh(test_bullet_rotated, style='wireframe', color='yellow', line_width=5)
pl2.show()


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
