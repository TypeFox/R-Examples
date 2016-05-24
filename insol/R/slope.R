slope = function(cgrad,degrees=FALSE){
    if (nargs() < 1) {
        cat("USAGE: slope(cellgradient, degrees=F) \n")
        return()
    }
	slope=acos(cgrad[,,3])
	if (degrees) slope=degrees(slope)
	return(slope)
}



