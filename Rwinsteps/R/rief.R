rief <- function(x, theta = seq(-4, 4, length = 100)){

  out <- riif(x, theta)
  out$e <- 1/sqrt(out$pq)

  class(out) <- "rief"

  return(out)
}
