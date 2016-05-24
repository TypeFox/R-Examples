lin.tailup <- function(dist.hydro, weight, parsil = parsil, range1 = range1)
{
	no <- length(dist.hydro[,1])
	np <- length(dist.hydro[1,])
	parsil*(matrix(rep(1, times = no*np), nrow = no) - dist.hydro/range1)*
		(dist.hydro < range1)*weight
}

