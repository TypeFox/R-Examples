getAutoGridSize <-
function(nL) {
	index10 = floor(log(nL, base=10))
	if(nL/(10^index10) < 3) {
		index10 = index10-1
	}
	grid.size = 10^(1:index10)
	return(grid.size)
}

