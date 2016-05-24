`predictMu` <-
function(gp, newTheta) {
	if (gp$constantMean == 0) {
		return (gp$Bhat%*%c(1, newTheta))
	}
	return (gp$mu[1])   ## otherwise mu is constant
}

