Cftheta <-
function(theta.s, lambda.r, times, delta, type.t, K, covar, m, s) {
  theta.s <- dnorm(theta.s, mean = 0, sd = 10) * exp(sum(covar[delta==1, s] * theta.s) - sum(lambda.r * m))
  return(theta.s)
}
