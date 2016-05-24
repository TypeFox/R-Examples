epa.tailup <- function(dist.hydro, weight, parsil = parsil, range1 = range1)
{
	no <- length(dist.hydro[,1])
	np <- length(dist.hydro[1,])
	parsil/(16*range1^5)*(dist.hydro - range1)^2*
		(16*range1^3*(dist.hydro*0 + 1) +
		17*range1^2*dist.hydro -
		2*range1*dist.hydro^2 - 
		dist.hydro^3)*
		(dist.hydro < range1)*weight
}

