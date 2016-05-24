`h.inv` <-
function(u1, u2, rho) {
  y <- pnorm(qnorm(u1)*sqrt(1-rho^2) + rho*qnorm(u2))
  return(y)
}

