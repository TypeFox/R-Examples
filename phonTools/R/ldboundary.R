# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

ldboundary = function (mean1, mean2, covariance, prior1 = .5, prior2 = .5, add = F, ...){
  w = solve (covariance) %*% (mean1-mean2)
  x0 = 1/2 * (mean1+mean2) - (log(prior1/prior2) /
       t(mean1-mean2) %*% solve(covariance) %*% (mean1-mean2)) %*%
       (mean1 - mean2)
  x0 = matrix (x0,2,1)
  coeffs = c(intercept = -(t(w) %*% x0)/-w[2],slope = -w[1]/w[2])
  if (add) abline (coeffs, ...)
  coeffs
}

