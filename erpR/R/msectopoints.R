	msectopoints<-function(a, lengthsegment, startmsec, endmsec){
	
	#lengthsegment=lengthsegment-1
	totmsec=endmsec-(startmsec) #total duration in msec
	msecstep=(lengthsegment-1)/totmsec #how many points is a msec. 
	# notice the - 1. This is because it doesn't make sense to say how many point is a msec.
	# what you actually are counting is how many points there are between two msec
	# each msec step is not single point, but the space between two points (a small segment between two msec). For this reason you have to subtract one.
	x=(a-startmsec)*msecstep
	return(x+1)
	# notice that you return x+1. This is because the points in your vector goes from 1 to length of the vector.
	# when you are at point 1. you are at step 0 (i.e., no "small segment"). The maximum number of segment you can reach is lengthsegment-1.
	
	}
