mar.taildown <- function(dist.hydro, a.mat, b.mat, parsil = parsil, 
	range1 = range1, useTailDownWeight, weight = NULL)
{
	np <- length(a.mat[1,])
	flow.connect <- b.mat == 0
	pos.h <- dist.hydro > 0
	a.eq.b <- a.mat == b.mat
	part1 <- 1*(flow.connect & !pos.h)
	part1[is.nan(part1)] <- 0
	part2 <- log(90*dist.hydro/range1 + 1)/(90*dist.hydro/range1)*(flow.connect & pos.h)
	part2[is.nan(part2)] <- 0
	part3 <- (log(90*a.mat/range1 + 1) - log(90*b.mat/range1 + 1))/((90*a.mat - 90*b.mat)/range1)*
			(!flow.connect & !a.eq.b)
	part3[is.nan(part3)] <- 0
	part4 <- (1/(90*a.mat/range1 + 1))*(!flow.connect & a.eq.b)
	part4[is.nan(part4)] <- 0
	V <- parsil*(part1 + part2 + part3 + part4)
	if(useTailDownWeight == TRUE) V <- V*weight
	V
}

