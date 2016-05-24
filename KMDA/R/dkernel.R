dkernel <-
function(x,y,rho){   
	dxy=mdist(x,y)
	out=exp(-dxy^2/rho)
	return(out)
}
