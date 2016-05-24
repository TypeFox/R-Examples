CUpdTheta <-
function(theta, lambda.r, times, delta, type.t, K, covar) {
  p <- length(theta)
  theta.r <- theta
  m <- CGaM(times, delta, type.t, K, covar, theta)
  for(s in 1:p){
    theta.str <- theta.r
    theta.str[s] <- rnorm(1, mean = theta[s], sd = sqrt(2)) 
    m.str <- CGaM(times, delta, type.t, K, covar, theta.str)
    a <- Cftheta(theta.str[s], lambda.r, times, delta, type.t, K, covar, m.str, s)
    b <- Cftheta(theta.r[s], lambda.r, times, delta, type.t, K, covar, m, s)
    pr <- min(1,  a / b)
    if(runif(1) <= pr) {
      theta[s] <- theta.str[s]
    }
  }
  return(theta)
}
