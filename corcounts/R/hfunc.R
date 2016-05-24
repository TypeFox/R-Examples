`hfunc` <-
function(u1, u2, rho) {
  y <- pnorm((qnorm(u1)-rho*qnorm(u2))/sqrt(1-rho^2))
  return(y)
}

