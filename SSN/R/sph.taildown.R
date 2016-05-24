sph.taildown <-
function(dist.hydro, a.mat, b.mat, parsil = parsil, range1 = range1,
	useTailDownWeight, weight = NULL)
{
	flow.connect <- b.mat == 0
	no <- length(a.mat[,1])
	np <- length(a.mat[1,])
	V <- parsil*(matrix(rep(1, times = no*np), nrow = no) - 1.5*(dist.hydro/range1) +
		0.5*(dist.hydro/range1)^3)*(a.mat < range1)*flow.connect +
		parsil*(matrix(rep(1, times = no*np), nrow = no) - 1.5*b.mat/range1 +
		0.5*a.mat/range1)*(matrix(rep(1, times = no*np), nrow = no) - a.mat/range1)^2*
		(a.mat < range1)*(1 - flow.connect)
	if(useTailDownWeight == TRUE) V <- V*weight
	V
}

