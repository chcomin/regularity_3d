''' Calculate the regularity of a 3D set of points using the method described 
	in "Quantifying the Regularity of a 3D Set of Points on the Surface of an 
	Ellipsoidal Object" (to be published). The main function to be called is
	regularity_3d().
'''

import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import minimize
from shapely.geometry import Polygon, Point
from igraph import Graph
import scipy.interpolate
import shelve

import util
import voronoi_3d
import hexagonal_grid
import visualization

def algebraic_ellipsoid_distance(pars, data):
	'''Calculate the algebraic distance for points in 'data' with respect
	   to the ellipsoid described by 'pars', which is a tuple containing
	   the center of the ellipsoid end the axes sizes (x0, y0, z0, a, b, c) '''

	x0, y0, z0, a, b, c = pars
	v = np.sum(np.abs((data[0]-x0)**2/a**2 + (data[1]-y0)**2/b**2 + (data[2]-z0)**2/c**2 - 1))
	return v	
	
def create_ellipsoid(x0, y0, z0, a, b, c):
	''' Create an ellipsoid for testing purposes and also for plotting functions '''

	theta = np.linspace(0, np.pi, 30)
	phi = np.linspace(0, 2*np.pi, 60)
	T, P = np.meshgrid(theta, phi)
	T, P = T.ravel(), P.ravel()
	r = 1./np.sqrt(np.sin(T)**2*np.cos(P)**2 + np.sin(T)**2*np.sin(P)**2 + np.cos(T)**2)
	X = r*np.sin(T)*np.cos(P)
	Y = r*np.sin(T)*np.sin(P)
	Z = r*np.cos(T)
	x_ellipsoid = a*X+x0
	y_ellipsoid = b*Y+y0
	z_ellipsoid = c*Z+z0
	points = np.array([x_ellipsoid, y_ellipsoid, z_ellipsoid]).T

	return points
	
def normalize_vec(vec):
	''' Create unit vector'''

	vec_magnitude = np.sqrt(np.sum(vec**2))	
	vec_normal = vec/vec_magnitude
	return vec_normal
	
def perpendicular(vet):
	'''Get versors that are orthogonal to vet'''
	
	c1_vet, c2_vet, c3_vet = vet
	
	if c1_vet==0:
		e1 = np.array([1, 0, 0])	
	else:
		c2_e1 = 1.
		c3_e1 = 0. 
		c1_e1 = -c2_vet*c2_e1/c1_vet
		e1 = np.array([c1_e1, c2_e1, c3_e1])
		
	e2 = np.cross(vet,e1)

	e1 = e1/np.sqrt(e1[0]**2+e1[1]**2+e1[2]**2)
	e2 = e2/np.sqrt(e2[0]**2+e2[1]**2+e2[2]**2)	
	
	return e1, e2	

def get_rel_angles(points, adjacency_list):
	''' Get angles between neighbors of each point in 'adjacency_list'.'''

	rel_angles = []
	for node_ref, neighbors in enumerate(adjacency_list):
		
		p_ref = points[node_ref]
		normal_plane = -normalize_vec(p_ref)
		
		e1, e2 = perpendicular(normal_plane)
		transf_mat = np.stack([e1, e2])
		
		angles = []
		for neighbor in neighbors:
			p_nei = points[neighbor]
			nei_vector = p_nei - p_ref	
			p_nei_projected = np.dot(transf_mat, nei_vector)
			angle_nei = np.arctan2(p_nei_projected[1], p_nei_projected[0])
			if angle_nei<0:
				angle_nei += 2*np.pi
			angles.append(angle_nei)
			
		k = np.argsort(angles)
		ordered_angles = np.array(angles)[k]
		ordered_neighbors = np.array(neighbors)[k]
		
		rel_angles.append([])
		for neighbor_index in range(len(ordered_neighbors)-1):
			rel_angle = ordered_angles[neighbor_index+1] - ordered_angles[neighbor_index]
			rel_angles[-1].append(rel_angle)
		rel_angle = (2*np.pi-ordered_angles[-1]) + ordered_angles[0]
		rel_angles[-1].append(rel_angle)
		
	return rel_angles
	
def order_neighbors(adjacency_list, points):
	'''Order neighbors of a node according to their spatial distance'''

	ordered_adj_list = []
	for node_ref, neighbors in enumerate(adjacency_list):
		
		p_ref = points[node_ref]
		normal_plane = -normalize_vec(p_ref)
		
		e1, e2 = perpendicular(normal_plane)
		transf_mat = np.stack([e1, e2])
		
		angles = []
		for neighbor in neighbors:
			p_nei = points[neighbor]
			nei_vector = p_nei - p_ref	
			p_nei_projected = np.dot(transf_mat, nei_vector)
			angle_nei = np.arctan2(p_nei_projected[1], p_nei_projected[0])
			if angle_nei<0:
				angle_nei += 2*np.pi
			angles.append(angle_nei)
			
		k = np.argsort(angles)
		ordered_neighbors = [neighbors[index] for index in k]
		# Set neighbor with smallest index as first item
		first_neighbor_index = np.argmin(ordered_neighbors)
		num_neighbors = len(ordered_neighbors)
		ordered_neighbors_smallest_first = []
		for index in range(num_neighbors):
			new_index = (first_neighbor_index+index)%num_neighbors
			ordered_neighbors_smallest_first.append(ordered_neighbors[new_index])
			
		ordered_adj_list.append(ordered_neighbors_smallest_first)
		
	return ordered_adj_list	
	
def polygonality(points, adjacency_list, region_border, ellipsoid_axes, ref_angle=np.pi/3):
	'''Calculate polygonality for a set of 3D points.

	   Parameters:
	   -----------
	   points : numpy array
	   		Nx3 array containing the position of N points.
	   adjacency_list : list of lists
	   		Adjacency list of a graph describing neighborhoods between the points
	   region_border : shapely Polygon
	   		Polygon describing the border of the point cloud
	   ellipsoid_axes : tuple
	   		Tuple containing the axes sizes of the ellipsoid (a, b, c)
	   ref_angle : float
	   		Reference angle to use for polygonality calculation
	'''

	avg_abs_angles_diff = get_ref_polygonality_dist(len(points), region_border, ellipsoid_axes, ref_angle)
	avg_grid_angles = np.mean(avg_abs_angles_diff)
	std_grid_angles = np.std(avg_abs_angles_diff)
	min_grid_angles = np.min(avg_abs_angles_diff)

	rel_angles = get_rel_angles(points, adjacency_list)
	N = len(rel_angles)
	avg_rel_angles = np.zeros(N)
	for i in range(N):
		rel_ang = rel_angles[i]
		avg_rel_angles[i] = np.mean(np.abs(np.array(rel_ang)-ref_angle))

	t = (avg_rel_angles-min_grid_angles)/std_grid_angles
	poly_ind = 1/(1+np.maximum(t, 0))

	return poly_ind, avg_rel_angles
	
def project_point_on_ellipsoid(point, a, b, c):
	''' Project a point onto the surface of an ellipsoid'''
		
	dist_r = np.sqrt(np.sum(point**2))
	versor_r = point/dist_r
	q = 1./np.sqrt(versor_r[0]**2/a**2 + versor_r[1]**2/b**2 + versor_r[2]**2/c**2)
	projected_point = q*versor_r

	return projected_point

def project_point_into_2D(point, a, b, c):
	''' Spherical projection of a point'''

	theta = np.arccos(point[2]/c)    
	phi = np.arctan2(a*point[1], b*point[0])

	return theta, phi

def project_2D_point_into_3D(theta, phi, a, b, c):
	''' Get Cartesian coordinates of points on the surface of an ellipsoid
		described by spherical coordinates. '''

	x = a*np.sin(theta)*np.cos(phi)
	y = b*np.sin(theta)*np.sin(phi)
	z = c*np.cos(theta)

	return x, y, z

def shift_2d_points(points_in_2d, ref_phi=None):
	''' Shift points so as to diminish the influence of spherical projection.
	    We want points to be far away from the sphere poles (small and large
		phi values).'''
	    

	if ref_phi is None:
		hist, bins = np.histogram(points_in_2d[:,1], np.linspace(-np.pi, np.pi, 100))
		ind = np.argmax(hist)
		ref_phi = bins[ind]
	
	points_in_2d_shifted = points_in_2d.copy()
	points_in_2d_shifted[:,1] = points_in_2d_shifted[:,1] - ref_phi
	ind = np.nonzero(points_in_2d_shifted[:,1]<=-np.pi)[0]
	points_in_2d_shifted[ind,1] = 2*np.pi + points_in_2d_shifted[ind,1]
	ind = np.nonzero(points_in_2d_shifted[:,1]>np.pi)[0]
	points_in_2d_shifted[ind,1] = points_in_2d_shifted[ind,1] - 2*np.pi

	return points_in_2d_shifted, ref_phi
	
def get_border_2D(points, ellipsoid_axes, shift_points=False):
	''' Obtain the border of the point cloud, described by array points,
	    along the surface of an ellipsoid (with axes set by ellipsoid_axes).'''

	a, b, c = ellipsoid_axes
	points_in_2d = np.zeros((len(points), 2))
	for p_index, point in enumerate(points):
		p_point = project_point_on_ellipsoid(point, a, b, c) 
		theta, phi = project_point_into_2D(p_point, a, b, c)
		points_in_2d[p_index] = (theta, phi)

	if shift_points:
		points_in_2d, ref_phi = shift_2d_points(points_in_2d)

	chull = ConvexHull(points_in_2d)
	chull_points = points_in_2d[chull.vertices]
	region_border = Polygon(chull_points)	

	if shift_points:
		return region_border, chull.vertices, ref_phi
	else:
		return region_border, chull.vertices

def get_border_3D(region_border, ellipsoid_axes):
	''' Project the border obtained by function get_border_2D()
		onto the surface of an ellipsoid'''

	a, b, c = ellipsoid_axes
	region_border_3d = []
	for p in np.linspace(0, 1, 100):
		pb = region_border.exterior.interpolate(p, True)
		x, y, z = project_2D_point_into_3D(pb.coords[0][0], pb.coords[0][1], a, b, c)
		region_border_3d.append((x, y, z))

	return region_border_3d
	
def get_best_triangle_size(num_points, region, ellipsoid_axes, size_limits=(8, 16)):
	''' Obtain the optimal number of points for generating the hexagonal grid inside
		the ellipsoidal region defined by 'region' and 'ellipsoid_axes'. We want a 
		parameter n for function hexagonal_grid.generate_grid(n) such that the generated
		grid has a number of points that is as close as possible to 'num_points'.'''

	a, b, c = ellipsoid_axes
	size_values = range(size_limits[0], size_limits[1]+1)
	smallest_diff = num_points
	for n in size_values:
		grid_pos, adjacency_list = hexagonal_grid.generate_grid(n=n)
		
		nodesInRegion = 0
		for i,p in enumerate(grid_pos):
			grid_p_ellip = project_point_on_ellipsoid(p, a, b, c)	
			theta, phi = project_point_into_2D(grid_p_ellip, a, b, c)
			if region.contains(Point(theta, phi)):
				nodesInRegion += 1

		if abs(nodesInRegion-num_points)<smallest_diff:
			smallest_diff = abs(nodesInRegion-num_points)
			best_n = n
			
	return best_n	
	
def get_grid_in_ellipsoid(grid_adjacency_list_all, grid_pos_3D, region_border, ellipsoid_axes):
	''' Project a grid of points, with adjacencies 'grid_adjacency_list_all' and positions
	    'grid_pos_3D' onto the surface of an ellipsoid delimited by 'region_border'.'''

	a, b, c = ellipsoid_axes
	edges = []
	for node1, neighbors in enumerate(grid_adjacency_list_all):
		for node2 in neighbors:
			if node1>node2:
				edges.append((node1, node2))
	g_grid_all = Graph(edges=edges)	

	grid_pos_ellip = []
	nodes2keep = []
	for i,p in enumerate(grid_pos_3D):
		grid_p_ellip = project_point_on_ellipsoid(p, a, b, c)	
		theta, phi = project_point_into_2D(grid_p_ellip, a, b, c)
		if region_border.contains(Point(theta, phi)):
			grid_pos_ellip.append(grid_p_ellip)
			nodes2keep.append(i)
	grid_pos_ellip = np.array(grid_pos_ellip)

	g_grid = g_grid_all.subgraph(nodes2keep)
	
	return g_grid, grid_pos_ellip
	
def get_ref_polygonality_dist(num_points, region_border, ellipsoid_axes, ref_angle=np.pi/3):
	''' Get reference angle distribution for a hexagonal grid. The average angle from 
		this distirbution is used in the definition of polygonality. '''

	a, b, c = ellipsoid_axes
	best_n = get_best_triangle_size(num_points, region_border, (a, b, c))
	grid_pos_3D, grid_adjacency_list_all = hexagonal_grid.generate_grid(n=best_n)
	g_grid, grid_pos_ellip = get_grid_in_ellipsoid(grid_adjacency_list_all, grid_pos_3D, region_border, (a, b, c))
	grid_adjacency_list = g_grid.get_adjlist()
			
	angles_grid = get_rel_angles(grid_pos_ellip, grid_adjacency_list)

	avg_abs_angles_diff = np.zeros(len(angles_grid))
	for point_index, l in enumerate(angles_grid):
		avg_abs_angles_diff[point_index] = np.mean(np.abs(np.array(l)-ref_angle))
			
	return avg_abs_angles_diff	

def regularity_3d(points, d_max=5., plot_data=False):
	'''Calculate the polygonality for a set of 3D points.

	   Parameters:
	   -----------
	   points : numpy array
	   		Nx3 array containing the position of N points.
	   d_max : float
	   		The maximum distance allowed in the Voronoi graph. Please refer to
			the paper for an explanation.
	   plot_data : bool
	   		If True, plots showing some steps of the method will be shown using
			the Plotly library
	'''

	# Fit ellipsoid to data, using the centroid as initial guess
	initial_x0, initial_y0, initial_z0 = np.mean(points, axis=0)
	initial_a, initial_b, initial_c = np.abs(np.max(points, axis=0) - np.min(points, axis=0))/2

	res = minimize(algebraic_ellipsoid_distance, 
			x0=(initial_x0, initial_y0, initial_z0, initial_a, initial_b, initial_c), args=(points.T))
	x0, y0, z0, a, b, c = res['x']

	points_trans = points - np.array([x0, y0, z0])

	# Create Voronoi graph
	g_points = voronoi_3d.create_voronoi_graph(points_trans, dist_thresh=d_max)
	adjacency_list = g_points.get_adjlist()

	# Get the border of the point cloud along the surface of the ellipsoid
	region_2D, border_vertices = get_border_2D(points_trans, (a, b, c))

	is_border = np.zeros(g_points.vcount(), dtype=np.uint8)
	is_border[border_vertices] = 1
	g_points.vs['is_border'] = is_border.tolist()

	# Calculate regularity of the graph
	poly, avg_rel_angles = polygonality(points_trans, adjacency_list, region_2D, (a, b, c), np.pi/3.)

	# Store the results 
	prop_dict = {'graph':g_points, 'points':points_trans,
				 'polygonality':poly, 'ellipsoid_center':(x0, y0, z0), 
				 'ellipsoid_axes':(a, b, c), 'border_polygon':region_2D}

	# Plot data, if desired
	if plot_data:
		visualization.generate_plots(g_points, points_trans, region_2D, (a,b,c), poly)

	return prop_dict

# An example of the application of the methodology is included below

if __name__=='__main__':

	input_file = 'DriedGnattData'
	output_file = 'measurements'

	facets = util.read_facets('data/%s_facet_pts.csv'%input_file)

	# Calculate facet centers
	facet_centers = util.get_facet_centers(facets)
	prop_dict = regularity_3d(facet_centers, d_max=5., plot_data=True)

	# Store the result in a database
	with shelve.open(output_file) as mea_db:
		mea_db[input_file] = prop_dict
