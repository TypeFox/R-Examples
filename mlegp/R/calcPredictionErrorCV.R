`calcPredictionErrorCV` <-
function(gp, X, newTheta, varMatrix, nugget) {
        r = calcCorOneObs(X, gp$beta, gp$a, newTheta)
        return (  (gp$sig2 + nugget) - gp$sig2*r%*%solve(varMatrix)%*%t(r)*gp$sig2  )
}

