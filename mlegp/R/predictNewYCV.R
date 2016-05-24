`predictNewYCV` <-
function(gp, X, Z, muMatrix, newTheta, varMatrix) {
  r = calcCorOneObs(X, gp$beta, gp$a, newTheta)
  newY = predictMu(gp, newTheta) + gp$sig2 * r %*% solve(varMatrix) %*% (Z - muMatrix)
  return (newY)
}

