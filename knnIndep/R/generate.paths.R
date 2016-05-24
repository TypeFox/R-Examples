generate.paths <-
function(index,rx,ry,N){
	
	x = rx - rx[index]
	y = ry - ry[index]
	x = pmin(abs(rx - rx[index]),N -abs(rx - rx[index]))
	y = pmin(abs(ry - ry[index]),N -abs(ry - ry[index]))
	#remove pt at origin
	x = x[-index]
	y = y[-index]
	
	dists = pmax(x,y)
	dsorted = sort(dists)
	return(dsorted)
}
