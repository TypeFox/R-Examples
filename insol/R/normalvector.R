normalvector <-
function(slope,aspect) {
	if (nargs() < 2 ) {cat("USAGE: normalvector(slope,aspect) \n"); return()}
	sloper = radians(slope)	
	aspectr = radians(aspect)
	nvx = sin(aspectr)*sin(sloper)
	nvy = -cos(aspectr)*sin(sloper)
	nvz = cos(sloper)
	return(cbind(nvx,nvy,nvz))
}

