exp.taildown <- function(dist.hydro, a.mat, b.mat, parsil = parsil, 
	range1 = range1, useTailDownWeight, weight = NULL)
{
	V <- parsil*exp(-3*dist.hydro/range1)
	if(useTailDownWeight == TRUE) V <- V*weight
	V
}

