posttrace <-
function(fmodel, y, zeta = seq(zmin, zmax, length = length), 
   zmin = -3, zmax = 3, length = 100, ...) 
{
	if (is.vector(y)) {
		y <- matrix(y, 1, length(y))
	}
	post <- rep(NA, length(zeta)) 
	for (i in 1:length(zeta)) {
		post[i] <- fmodel(zeta[i], y, ...)$post	
	}
	return(list(zeta = zeta, post = post))	
}