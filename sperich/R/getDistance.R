getDistance <-
function(point.a, point.b, resolution=1){
	distance.ab <- 0
	
	#get distances (grid)
	dx <- point.a[1]-point.b[1]
	dy <- point.a[2]-point.b[2]

	#get distances (degree)
	xdist <- dx * resolution
	ydist <- dy * resolution

	#pythagoras
	distance.ab <- sqrt(xdist^2 + ydist^2)

	return(distance.ab)
}