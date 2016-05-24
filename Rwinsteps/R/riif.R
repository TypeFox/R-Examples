riif <- function(x, theta = seq(-4, 4, length = 100)){

  out <- rirf(x, theta)
  out$pq <- out$p*(1 - out$p)

  class(out) <- "riif"

  return(out)
}
