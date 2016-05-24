pointstomsec=function(a, lengthsegment, startmsec, endmsec)
{
	totmsec=endmsec-(startmsec) #total duration in msec
	pointsstep=totmsec/(lengthsegment-1) #how many msec is a point.
	# notice the - 1. This is because it doesn't make sense to say how many msec there are in a point.
	# what you actually are counting is how many msec there are between two points.
	# each point step is the space between two points (a small segment between two points). For this reason you have to subtract one.	
	x=((a-1)*pointsstep)+(startmsec)
	# notice that you check a-1. This is because the points in your vector goes from 1 to length of the vector. 
	# when you are at point 1. you are at step 0 (i.e., no "small segment"). The maximum number of segment you can reach is lengthsegment-1.

	return(x)

}
