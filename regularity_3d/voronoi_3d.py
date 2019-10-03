''' Generation of a Delaunay graph with some additional contraints, as described in "Quantifying 
	the Regularity of a 3D Set of Points on the Surface of an Ellipsoidal Object" (to 
	be published).'''

import numpy as np
from igraph import Graph
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, Point
import util

def create_ellipsoid(x0, y0, z0, a, b, c):
	''' Create an ellipsoid for testing purposes'''

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
	
def create_voronoi_graph(points, dist_thresh, infty_factor=3):
	''' Create Delaunay graph considering additional constraints on the allowed
	    connections (thus the name Voronoi graph). Parameter dist_thresh sets the 
		maximum distance allowed (d_{max} in the article). Parameter infty_factor
		is used for creating points far away from the cloud of points in order to
		have a closed Voronoi cell for every point in 'points'. It has no influence
		in the result provided it is large enough. '''

	min_pos = np.min(points, axis=0)
	max_pos = np.max(points, axis=0)
	L = np.max(max_pos-min_pos)
	max_v = np.max(max_pos)
	min_v = np.max(min_pos)
	max_infty = max_v+infty_factor*L
	min_infty = min_v-infty_factor*L
	points_at_infinity = np.array([(min_infty, min_infty, min_infty), (min_infty, max_infty, min_infty),
								   (max_infty, max_infty, min_infty), (max_infty, min_infty, min_infty),
								   (min_infty, min_infty, max_infty), (min_infty, max_infty, max_infty),
								   (max_infty, max_infty, max_infty), (max_infty, min_infty, max_infty)])
	points_with_infinity = np.concatenate((points, points_at_infinity), axis=0)
	index_infty_points = set(range(len(points), len(points)+8))

	vor = Voronoi(points_with_infinity)

	edges = []
	distances = []
	for i in range(len(vor.ridge_points)):
		nodes = vor.ridge_points[i]
		vertices = vor.ridge_vertices[i]
		if nodes[0] not in index_infty_points and nodes[1] not in index_infty_points:
			# We want to get the smallest distance between the middle point in 'nodes' and 
			# their respective Voronoi ridge
			# First, define new coordinate on ridge plane. Origin of the ridge plane is the middle 
			# point between original points
			pos_node1 = vor.points[nodes[0]]
			pos_node2 = vor.points[nodes[1]]
			pos_middle = (pos_node1+pos_node2)/2.	
			# First axis is the ridge normal
			ridge_normal = pos_node2-pos_middle
			ridge_normal = ridge_normal/np.sqrt(np.sum(ridge_normal**2))
			# Second axis goes from middle point to some vertex on the ridge border
			v1 = vor.vertices[vertices[0]] - pos_middle
			v1 = v1/np.sqrt(np.sum(v1**2))
			# Third axis is the cross product between normal and v1
			v2 = np.cross(ridge_normal, v1)
			M = np.array([v1, v2, ridge_normal]).T
			Mi = np.linalg.inv(M)
			# Get coordinates of ridge vertices on the ridge plane coords
			ridge_coords = np.dot(Mi, (vor.vertices[vertices] - pos_middle).T).T
			poly = Polygon(ridge_coords[:,:2])
			dist_middle_nodes = np.sqrt(np.sum((pos_node2-pos_middle)**2))
			if poly.contains(Point(0, 0)):
				dist_middle_ridge = 0.
			else:
				dist_middle_ridge = poly.exterior.distance(Point(0, 0))
				
			if dist_middle_ridge<=dist_thresh:
				edges.append(list(nodes))
				distances.append(dist_middle_ridge)					
				
				
	g = Graph(edges=edges)
	
	return g
	
if __name__=='__main__':

	# Test the method on artificial data or on the real data	
	#points = create_ellipsoid(0, 0, 0, 3, 1, 1)	
	#X, Z, Y = np.meshgrid([0, 1, 2], [0, 1, 2], [0, 1, 2])
	#points = np.array([X.ravel(), Y.ravel(), Z.ravel()]).T

	facets = util.read_facets('DriedGnattFemaleData.csv')
	facet_centers = []
	for facet in facets:
		facet_center = np.mean(facet, axis=0)
		facet_centers.append(facet_center)
	facet_centers = np.array(facet_centers)
	points = facet_centers.copy()

	g = create_voronoi_graph(points, dist_thresh=5., infty_factor=3)

			
	