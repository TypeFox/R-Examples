postmode <-
function(fmodel, y, zmin = -5, zmax = 5, ...) 
{
	tmp <- optimize(function(zeta, ...) fmodel(zeta, ...)$post, 
		interval = c(zmin, zmax), maximum = TRUE, y = y, ...)
	names(tmp) <- c("zeta", "post")
	return(tmp)
}