alpha.divergence <- function(dist1, dist2)
{	
	warning("Use of alpha.divergence(...) is deprecated!")
	return(alphaDivergence(dist1,dist2))
}

alphaDivergence <-
function(dist1, dist2)
{
	return ( 1-sum(sqrt(dist1*dist2)) )
	
}
