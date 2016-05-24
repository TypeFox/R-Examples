

hellinger.distance <-
function(x, y)
{
	warning("Use of hellinger.distance(...) is deprecated!")
	return(hellingerDistance(x,y))
}

hellingerDistance <-
function(x, y)
{
	return( 1/(sqrt(2))* sqrt( sum( (sqrt(x)-sqrt(y))**2 ) ))
	
}
