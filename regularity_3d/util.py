''' Some utility functions specific to the fundus gnat data '''

import numpy as np
import csv

def read_facets(filename):
    '''Read facets from a given .csv file'''

    data = csv.reader(open(filename, 'r'))
    line = data.__next__()
    curr_index = int(line[0])
    pos = list(map(float, line[1:]))
    facets = []
    facet = [pos]
    for line in data:
        facet_index = int(line[0])
        pos = list(map(float, line[1:]))

        if facet_index==curr_index:
            facet.append(pos)
        else:
            facets.append(np.array(facet))
            facet = [pos]
            curr_index = facet_index
    facets.append(np.array(facet)) 
    return facets
    
def get_facet_centers(facets):
	''' Calculate the centroid of the fundus gnat eye facets'''

	facet_centers = []
	for facet in facets:
		facet_center = np.mean(facet, axis=0)
		facet_centers.append(facet_center)
	facet_centers = np.array(facet_centers)
	return facet_centers

def plot_all_facets(facets, facet_colors):
	''' Plot facets of the fundus gnat data''' 

	x = []
	y = []
	z = []
	c = []
	for facet_index, facet in enumerate(facets):
		x.extend(facet[:,0])
		y.extend(facet[:,1])
		z.extend(facet[:,2])
		c.extend([facet_colors[facet_index]]*len(facet))

	layout = go.Layout(
		margin=dict(
			l=0,
			r=0,
			b=0,
			t=0
		),
	)
		
	facets = go.Scatter3d(
		x=x,
		y=y,
		z=z,
		mode='markers',
		marker=dict(
			size=1,
			color=c,                # set color to an array/list of desired values
			opacity=1.,
			colorscale='Viridis',
			colorbar=dict(
				title='Regularity'
			),

		),
		showlegend=False
	)		

	data = [facets] 
	fig = go.Figure(data=data, layout=layout)
	fig['layout'].update(height=1000, width=1000)
	plotly.offline.plot(fig, filename='full_facets.html')
	