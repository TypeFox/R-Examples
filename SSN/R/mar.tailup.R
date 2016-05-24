mar.tailup <- function(dist.hydro, weight, parsil = parsil, range1 = range1)
{
	no <- length(dist.hydro[,1])
	np <- length(dist.hydro[1,])
	CovMat <- matrix(0, ncol = np, nrow = no)
	ind <- (1e-9 < dist.hydro/range1)
	CovMat[ind] <- parsil*log(90*dist.hydro[ind]/range1 + 1)/(90*dist.hydro[ind]/range1)
	ind <- 1e-9 >= dist.hydro/range1
	CovMat[ind] <-parsil
	CovMat*weight
}

