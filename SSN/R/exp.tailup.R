exp.tailup <- function(dist.hydro, weight, parsil = parsil, range1 = range1)
{
	parsil*exp(-3*dist.hydro/range1)*weight
}

