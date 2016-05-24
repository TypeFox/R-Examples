`calcPredictionError` <-
function(gp, newData, nugget = 0) {
        r = calcCorOneObs(gp$X, gp$beta, gp$a, newData)
        v =  gp$sig2 + nugget - gp$sig2*r%*%gp$invVarMatrix%*%t(r)*gp$sig2 
	if (v < 0) return (0)
	return (sqrt(v)) 
}

