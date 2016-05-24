rtrf <- function(x, theta = seq(-4, 4, length = 100)){

  out <- list(theta = theta)
  out$p <- apply(rbind(rirf(x, theta)$p), 1, sum)

  class(out) <- "rtrf"

  return(out)
}
