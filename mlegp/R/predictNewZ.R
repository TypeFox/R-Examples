`predictNewZ` <-
function(gp, newTheta) {
  r = calcCorOneObs(gp$X, gp$beta, gp$a, newTheta)
  newY = predictMu(gp, newTheta) + gp$sig2 * r %*% gp$invVarMatrix %*% (gp$Z - gp$mu)
  return (newY)
}

