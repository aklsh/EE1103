# Downsampling and Interpolation

    Create a table of values for x and y=tan(x), for x = 1, 1.1, 1.2 and 1.3
    Use this table of values to interpolate and find the value for y at x=1.15 and x=1.35
    Compare your interpolated values using different methods of interpolation and against  the true value of the function.


The attached file has y-data separated in steps of 25ps. Downsample it to 250 ps (take every 10th point) and plot it in gnuplot. Now, generate the intermediate points  using interpolation, and compare the root mean square error between the actual data set and the interpolated points.

 
	
out1_test0.csv
