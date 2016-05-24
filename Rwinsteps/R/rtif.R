rtif <- function(x, theta = seq(-4, 4, length = 100)){

  out <- list(theta = theta)
  out$pq <- apply(rbind(riif(x, theta)$pq), 1, sum)

  class(out) <- "rtif"

  return(out)
}
