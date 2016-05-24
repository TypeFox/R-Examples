"epi.dgamma" <- function(rr, quantiles = c(0.05, 0.95)){
   fold.variation <- rr[2]/rr[1]
   low.p <- abs(qnorm(quantiles[1], mean = 0, sd = 1))
   up.p <- abs(qnorm(quantiles[2], mean = 0, sd = 1))
   p <- low.p + up.p
   tau <- (p^2) / (log(fold.variation) * log(fold.variation))
   return(tau)
}
