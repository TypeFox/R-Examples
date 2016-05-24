

squared.hellinger.distance <-
function(x, y)
{
	warning("Use of squared.hellinger.distance(...) is deprecated!")
	return(squaredHellingerDistance(x,y))
}

squaredHellingerDistance <-
function(x, y)
{
       return (0.5*sum( (sqrt(x)-sqrt(y))**2  ))
}
