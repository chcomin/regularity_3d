''' Generation of an hexagonal grid on the surface of an ellipsoid '''

from math import sin, cos, sqrt, pi, atan2
import numpy as np
from scipy.spatial.distance import pdist, squareform, cdist

def barycentricCoords(p):
	''' Barycentric coords for triangle (-0.5,0),(0.5,0),(0,sqrt(3)/2)'''

	x,y = p
	l3 = y*2./sqrt(3.)
	l2 = x + 0.5*(1 - l3)
	l1 = 1 - l2 - l3
	return l1,l2,l3

def scalProd(p1,p2):
	return sum([p1[i]*p2[i] for i in range(len(p1))])

def slerp(p0,p1,t):
	'''# Uniform interpolation of arc defined by p0, p1 (around origin)
	   # t=0 -> p0, t=1 -> p1'''

	assert abs(scalProd(p0,p0) - scalProd(p1,p1)) < 1e-7
	ang0Cos = scalProd(p0,p1)/scalProd(p0,p0)
	ang0Sin = sqrt(1 - ang0Cos*ang0Cos)
	ang0 = atan2(ang0Sin,ang0Cos)
	l0 = sin((1-t)*ang0)
	l1 = sin(t    *ang0)
	return tuple([(l0*p0[i] + l1*p1[i])/ang0Sin for i in range(len(p0))])

def mapGridpoint2Sphere(p,s1,s2,s3):
	''' Map 2D point p to spherical triangle s1,s2,s3 (3D vectors of equal length)'''

	l1,l2,l3 = barycentricCoords(p)
	if abs(l3-1) < 1e-10: return s3
	l2s = l2/(l1+l2)
	p12 = slerp(s1,s2,l2s)
	return slerp(p12,s3,l3)
	
def generate_grid(n):
	''' Generate a Geodesic polyhedron. Parameter n sets the number of vertices 
	    at each side of the initial equilateral triangle. Thus, it sets the number
		of points in the final grid.'''

	theta = np.pi/3.
	phi = 3*np.pi/10.
	triangle_corners = np.array([(-0.5,0),(0.5,0),(0,sqrt(3)/2)])
	triangle_side = triangle_corners[1,0] - triangle_corners[0,0]
	
	L = triangle_side/((n-1)*np.tan(theta))
	hf = L*np.tan(theta)/2.
	hv = L/(2*np.cos(theta))
	hp = L*np.cos(theta)

	y = triangle_corners[0,1]
	first_x = triangle_corners[0,0]
	hexa_triang = []
	hexa_triang_parts = {'base':[], 'left':[], 'right':[], 'middle':[]}
	for i in range(n-1):
		num_points = n-i
		x = np.linspace(first_x, first_x+(num_points-1)*2*hf, num_points)
		hexa_triang.extend(list(zip(x, [y]*num_points)))
		
		if i==0:
			hexa_triang_parts['base'] = list(zip(x[1:-1], [y]*(num_points-2)))
		else:
			hexa_triang_parts['left'].append((x[0], y))
			hexa_triang_parts['right'].append((x[-1], y))
			hexa_triang_parts['middle'].extend(list(zip(x[1:-1], [y]*(num_points-2))))

		y += L+hp
		first_x += hf

	hexa_triang.append((first_x, y-L/2-hp+hv))
	hexa_triang = np.array(hexa_triang)
		
	s,c = 2/sqrt(5),1/sqrt(5)
	topPoints = [(0,0,1)] + [(s*cos(i*2*pi/5.), s*sin(i*2*pi/5.), c) for i in range(5)]
	bottomPoints = [(-x,y,-z) for (x,y,z) in topPoints]
	icoPoints = np.array(topPoints + bottomPoints)
	icoTriangs = [(0,i+1,(i+1)%5+1) for i in range(5)] +\
				 [(6,i+7,(i+1)%5+7) for i in range(5)] +\
				 [(i+1,(i+1)%5+1,(7-i)%5+7) for i in range(5)] +\
				 [(i+1,(7-i)%5+7,(8-i)%5+7) for i in range(5)]  

	grid_pos = []			 
	for triang in icoTriangs:
		s1, s2, s3 = icoPoints[list(triang)]
		
		for p in hexa_triang:
			new_point = mapGridpoint2Sphere(p, s1, s2, s3)
			grid_pos.append(new_point)
	grid_pos = np.array(grid_pos)	
	for k,v in hexa_triang_parts.items():
		hexa_triang_parts[k] = np.array(v)
		
	D = squareform(pdist(grid_pos))
	points_to_remove = set()
	for i in range(D.shape[0]):
		ind = np.nonzero(D[i,i+1:]<1e-8)[0]
		for j in ind:
			points_to_remove.add(i+1+j)
	points_to_keep = set(range(D.shape[0])) - points_to_remove
	final_grid_pos = grid_pos[np.array(list(points_to_keep))]
	
	D_ico = cdist(icoPoints, final_grid_pos)
	icoVertices = np.argmin(D_ico, axis=1)
	D = squareform(pdist(final_grid_pos))
	adjacency_list = []
	for i in range(D.shape[0]):
		k = np.argsort(D[i])
		if i in icoVertices:
			degree = 5
		else:
			degree = 6
		neighbors = k[1:degree+1]
		adjacency_list.append(neighbors.tolist())
		
	return final_grid_pos, adjacency_list

def generate_dual_grid(n):
	''' Generate a Goldberg polyhedron. Parameter n sets the number of vertices 
	    at each side of the initial equilateral triangle. Thus, it sets the number
		of points in the final grid.'''

	s = 1./n
	theta = np.pi/3.
	L = s/(2*np.sin(theta))
	h = L*np.cos(theta)
	h_pent = L/(2*np.cos(3*np.pi/10.))	
	triangle_corners = np.array([(-0.5,0),(0.5,0),(0,sqrt(3)/2)])

	y = triangle_corners[0,1] + L/2.

	first_x = triangle_corners[0,0] + s/2.
	hexa_triang = []
	for i in range(n-1):
		num_points = n-i
		x = np.linspace(first_x, first_x+(num_points-1)*s, num_points)
		hexa_triang.extend(list(zip(x, [y]*num_points)))
	 
		y += h
		first_x += s/2.
		x = np.linspace(first_x, first_x+(num_points-2)*s, num_points-1)
		hexa_triang.extend(list(zip(x, [y]*(num_points-1))))
		 
		y += L
	hexa_triang.append((first_x, y))
	hexa_triang = np.array(hexa_triang)
		
	s,c = 2/sqrt(5),1/sqrt(5)
	topPoints = [(0,0,1)] + [(s*cos(i*2*pi/5.), s*sin(i*2*pi/5.), c) for i in range(5)]
	bottomPoints = [(-x,y,-z) for (x,y,z) in topPoints]
	icoPoints = np.array(topPoints + bottomPoints)
	icoTriangs = [(0,i+1,(i+1)%5+1) for i in range(5)] +\
				 [(6,i+7,(i+1)%5+7) for i in range(5)] +\
				 [(i+1,(i+1)%5+1,(7-i)%5+7) for i in range(5)] +\
				 [(i+1,(7-i)%5+7,(8-i)%5+7) for i in range(5)]  

	grid_pos = []			 
	for triang in icoTriangs:
		s1, s2, s3 = icoPoints[list(triang)]
		
		for p in hexa_triang:
			new_point = mapGridpoint2Sphere(p, s1, s2, s3)
			grid_pos.append(new_point)
	grid_pos = np.array(grid_pos)	
	
	return grid_pos, hexa_triang
	
def project_point_on_ellipsoid(point, a, b, c):
	'''Project a point onto the surface of an ellipsoid'''
		
	dist_r = np.sqrt(np.sum(point**2))
	versor_r = point/dist_r
	q = 1./np.sqrt(versor_r[0]**2/a**2 + versor_r[1]**2/b**2 + versor_r[2]**2/c**2)
	projected_point = q*versor_r		
	return projected_point	

if __name__=='__main__':

	n = 11
	grid_pos, adjacency_list = generate_grid(n)	
		
	grid_pos_proj = np.zeros_like(grid_pos)
	for i,p in enumerate(grid_pos):
		grid_pos_proj[i] = project_point_on_ellipsoid(p, 1.1, 0.9, 0.9)	

	
