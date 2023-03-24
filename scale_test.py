#%%
"""
Testing the scaling functionality of pyvista
Author: Joseph Ko
Created: 3/23/23
Sources:
https://docs.pyvista.org/api/core/_autosummary/pyvista.DataSetFilters.transform.html#transform
"""
import pyvista as pv
import numpy as np

#%%
# create a rosette bullet

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
pl.add_mesh(bullet, show_edges=None, color = obj_color, opacity=op)
pl.show()

#%%
# randomly modify aspect ratio
sf_basal = 1.25 # basal scaling factor
sf_prism = 1.5 # prism scaling factor
t_mat = np.array([[sf_basal, 0, 0, 0],
                  [0,sf_basal, 0, 0],
                  [0, 0, sf_prism, 0],
                  [0, 0, 0, 1]])
bullet_modified = bullet.transform(t_mat, inplace=False)
pl.add_mesh(bullet_modified, style='wireframe', color='yellow', line_width=5)


#%%
copy_test = bullet.copy()
copy_test.plot()
# %%
print(bullet)
print(copy_test)
# %%
