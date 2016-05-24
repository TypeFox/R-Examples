sph.Euclid <- function(parsil, distance.matrix)
{
	CovMat <- (1 - 1.5*distance.matrix + 0.5*distance.matrix^3)
	CovMat[distance.matrix > 1] <- 0
	parsil*CovMat
}

