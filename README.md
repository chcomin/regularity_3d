# Regularity 3D

Calculate the regularity of a point cloud using the methodology described in "Quantifying the Regularity of a 3D Set of Points on the Surface of an Ellipsoidal Object" (to be published). It is assumed that the points are distributed approximately along the surface of an ellipsoid.

The main function to use is regularity_3d() from file [regularity_3d](regularity_3d/regularity_3d.py). Please refer to the documentation of that function (and also of the function polygonality() in the same file) for usage details. 

An example dataset of the fundus gnat eye is included in the data directory.

Example usage:

```python
  import regularity_3d
  
	input_file = 'DriedGnattData'
	output_file = 'measurements'

	facets = util.read_facets('data/%s_facet_pts.csv'%input_file)

	# Calculate facet centers
	facet_centers = util.get_facet_centers(facets)
	prop_dict = regularity_3d.regularity_3d(facet_centers, d_max=5., plot_data=True)

	# Store the result in a database
	with shelve.open(output_file) as mea_db:
		mea_db[input_file] = prop_dict
```


### Dependencies (version)
* Python version 3.7.1
* numpy (1.16.2)
* scipy (1.2.1)
* python-igraph (0.7.1.post6)
* shapely (1.6.4.post1)
* plotly (3.7.0)
