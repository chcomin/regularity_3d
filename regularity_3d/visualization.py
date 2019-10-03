''' This file contains some auxiliary plotting functions using 
	the Plotly library.
'''

import plotly
import plotly.graph_objs as go
import numpy as np
import hexagonal_grid
import regularity_3d

plotly.tools.set_credentials_file(username='chcomin', api_key='7jfNCE45TY6tSBwIXb5A')

def plot_data(lattice_p, points_p, el_points, border_points, filename):
	''' Plot on the same axes three point clouds. The visual properties of the points consider
		that the variables are respective to:
		lattice_p: a regular grid of points
		points_p: a set of points scattered along the surface of an ellipsoid
		el_points: points representing an ellipsoid
		In addition, the border of the region is represented as a set of points
		contained in 'border_points'. The plot is saved in file 'filename' '''
		
	x_lattice, y_lattice, z_lattice = lattice_p
	x_points, y_points, z_points = points_p
	x_ellipsoid, y_ellipsoid, z_ellipsoid = el_points
	x_border, y_border, z_border = border_points

	layout = go.Layout(
		margin=dict(
			l=0,
			r=0,
			b=0,
			t=0
		),
	)

	points = go.Scatter3d(
		x=x_points,
		y=y_points,
		z=z_points,
		mode='markers',
		marker=dict(
			size=4,
			color='red',                # set color to an array/list of desired values
			opacity=1.
		),
		showlegend=False
	)	
	
	lattice = go.Scatter3d(
		x=x_lattice,
		y=y_lattice,
		z=z_lattice,
		mode='markers',
		marker=dict(
			size=7,
			color='blue',                # set color to an array/list of desired values
			opacity=1.
		),
		showlegend=False,
	)

	ellipsoid = go.Scatter3d(
		x=x_ellipsoid,
		y=y_ellipsoid,
		z=z_ellipsoid,
		mode='markers',
		marker=dict(
			size=2,
			color='gray',                # set color to an array/list of desired values
			opacity=1.
		),
		showlegend=False
	)
	
	border = go.Scatter3d(
		x=x_border,
		y=y_border,
		z=z_border,
		mode='markers',
		marker=dict(
			size=2,
			color='green',                # set color to an array/list of desired values
			opacity=1.
		),
		showlegend=False
	)		
	data = [lattice, points, ellipsoid, border] 
	fig = go.Figure(data=data, layout=layout)
	plotly.offline.plot(fig, filename=filename)		
	
def plot_graph(g, points, filename, node_color='blue', cmap='Viridis', plot=True, cmin=None, cmax=None):
	''' Plot a graph g (an igraph data structure) having vertices with 
	    position given by the array 'points' '''

	Xn, Yn, Zn = points.T
	Xe=[]
	Ye=[]
	Ze=[]
	for e in g.get_edgelist():
		Xe += [Xn[e[0]],Xn[e[1]], None]
		Ye += [Yn[e[0]],Yn[e[1]], None]
		Ze += [Zn[e[0]],Zn[e[1]], None]
	
	layout = go.Layout(
		margin=dict(
			l=0,
			r=0,
			b=0,
			t=0
		),
	)

	edges = go.Scatter3d(
		x=Xe,
		y=Ye,
		z=Ze,
		mode='lines',
		line=dict(
			color='blue',       # set color to an array/list of desired values
		),
		showlegend=False
	)

	marker_dict = dict(
			size=6,
			color=node_color,          # set color to an array/list of desired values
			opacity=1.,
	)	
	if not isinstance(node_color, str):
		marker_dict['colorbar'] = dict(title='Polygonality')
		marker_dict['colorscale'] = cmap

		if cmax is None:
			marker_dict['cmax'] = max(node_color)
		else:
			marker_dict['cmax'] = cmax
		if cmin is None:
			marker_dict['cmin'] = min(node_color)
		else:
			marker_dict['cmin'] = cmin
			
	nodes = go.Scatter3d(
		x=Xn,
		y=Yn,
		z=Zn,
		mode='markers',
		marker=marker_dict,
		showlegend=False
	)
	
	if plot:
		data = [nodes, edges] 
		fig = go.Figure(data=data, layout=layout)
		fig['layout'].update(height=800, width=800)
		plotly.offline.plot(fig, filename=filename)			
	return nodes, edges
	
def generate_plots(graph, points, region_border, ellipsoid_axes, regularity, filename_prefix='plot'):
	''' Generate plots showing the results of some steps of the methodology. graph is an igraph graph
	    describing the adjacencies between the points. points is an Nx3 array containing the points to
		be plotted. region_border is a shapely polygon containing the border of the point cloud. 
		ellipsoid_axes is a tuple containing the sizes of the axes. regularity sets the colors of the
		points.'''

	a, b, c = ellipsoid_axes
	n = regularity_3d.get_best_triangle_size(len(points), region_border, (a, b, c))
	# Get grid with similar number of points as the array 'points'
	grid_pos_3D, grid_adjacency_list_all = hexagonal_grid.generate_grid(n=n)
	g_grid, grid_pos_ellip = regularity_3d.get_grid_in_ellipsoid(grid_adjacency_list_all, grid_pos_3D, region_border, (a, b, c))
	grid_adjacency_list = g_grid.get_adjlist()

	el_points = regularity_3d.create_ellipsoid(0, 0, 0, a, b, c)
	region_border_3d = regularity_3d.get_border_3D(region_border, (a, b, c))
	filename = '%s_points.html'%filename_prefix
	plot_data(grid_pos_ellip.T, points.T, el_points.T, 
			np.array(region_border_3d).T, filename)

	filename = '%s_grid_graph.html'%filename_prefix
	nodes1, edges1 = plot_graph(g_grid, grid_pos_ellip, filename, 
						plot=False)
	filename = '%s_voronoi_graph.html'%filename_prefix
	nodes2, edges2 = plot_graph(graph, points, filename, plot=False)
	fig = plotly.tools.make_subplots(rows=1, cols=2, specs=[[{'is_3d': True}, {'is_3d': True}]])
	fig.append_trace(nodes1, 1, 1)
	fig.append_trace(edges1, 1, 1)
	fig.append_trace(nodes2, 1, 2)
	fig.append_trace(edges2, 1, 2)

	fig['layout'].update(height=800, width=1400)
	filename = '%s_grid_and_voronoi_graph.html'%filename_prefix
	plotly.offline.plot(fig, filename=filename)	
	
	filename = '%s_points_regularity.html'%filename_prefix
	_, _ = plot_graph(graph, points, filename, 
					node_color=regularity, cmap='Viridis', plot=True)