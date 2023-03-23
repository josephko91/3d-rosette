#%%
# Function to return randomly sampled points on a spherical cap of a 2-sphere (i.e. 3-d sphere)
# Author: Joseph Ko
# Created: 3/16/23
# Sources: 
# https://stackoverflow.com/questions/38997302/create-random-unit-vector-inside-a-defined-conical-region/39003745#39003745
# https://math.stackexchange.com/questions/56784/generate-a-random-direction-within-a-cone/205589#205589
import pyvista as pv
import numpy as np

def random_spherical_cap(cone_angle_deg, cone_direction, num_points):
    cone_angle_rad = cone_angle_deg*(np.pi/180)
    # generate points on spherical cap centered at north pole
    z = np.random.uniform(np.cos(cone_angle_rad), 1, num_points)
    phi = np.random.uniform(0, 2*np.pi, num_points)
    x = np.sqrt(1-z**2)*np.cos(phi)
    y = np.sqrt(1-z**2)*np.sin(phi)
    points = np.column_stack((x, y, z))
    return points

cap_points = random_spherical_cap(30, 0, 50)
# print(cap_points)
sphere = pv.Sphere(radius=1, center=(0, 0, 0), direction=(0, 0, 1), 
                   theta_resolution=30, phi_resolution=30, start_theta=0, 
                   end_theta=360, start_phi=0, end_phi=180)
pl = pv.Plotter()
pl.background_color = 'black'
pl.add_mesh(sphere, show_edges=None, color = 'white', opacity=0.9)
pl.add_points(cap_points, color='red')

# rotate points
def norm_rows(v):
    if np.all(v==0):
        v_unit = np.array([1,0,0])
    else:
        if v.ndim == 1:
            v_norm = np.linalg.norm(v)
            v_unit = v/v_norm
        
        else:
            v_norm = np.linalg.norm(v, axis=1)
            v_unit = v/v_norm[:,None]
    return v_unit

north_vector = np.array([0, 0, 1])
cone_direction = np.array([-3, -4, 2])
cone_direction_norm = norm_rows(cone_direction)
u = norm_rows(np.cross(north_vector, cone_direction_norm)) # rotation axis
rot = np.arccos(np.dot(cone_direction_norm, north_vector)) # rotation angle in radians
ux = u[0]
uy = u[1]
uz = u[2]
# define rotation matrix
r11 = np.cos(rot) + (ux**2)*(1 - np.cos(rot))
r12 = ux*uy*(1 - np.cos(rot)) - uz*np.sin(rot)
r13 = ux*uz*(1 - np.cos(rot)) + uy*np.sin(rot)
r21 = uy*ux*(1 - np.cos(rot)) + uz*np.sin(rot)
r22 = np.cos(rot) + (uy**2)*(1 - np.cos(rot))
r23 = uy*uz*(1 - np.cos(rot)) - ux*np.sin(rot)
r31 = uz*ux*(1 - np.cos(rot)) - uy*np.sin(rot)
r32 = uz*uy*(1 - np.cos(rot)) + ux*np.sin(rot)
r33 = np.cos(rot) + (uz**2)*(1 - np.cos(rot))
rot_mat = np.array([[r11, r12, r13], 
                    [r21, r22, r23], 
                    [r31, r32, r33]])
# print(u)
# print(rot)
# print(rot_mat)
# print(cap_points)
cap_points_rot = np.matmul(rot_mat, cap_points.T)
cap_points_rot = cap_points_rot.T
pl.add_points(cap_points_rot, color='blue')
pl.show()


# %%
