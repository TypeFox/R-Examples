

calcModC01 <- function(X) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Compute and rank the 0-1 correlation
  # INPUT
  #   X: n x p matrix
  # OUTPUT
  #   phi: p*(p-1)/2 x 3 array
  # CALLED BY
  #  calcLarsNEW
  # ---------------------------------------------------------------
  p <- dim(X)[2]
  n <- dim(X)[1]
  rho <- cor(X)
  phi <- array(0,c(p*(p-1)/2,3))
  res <- .C("GGMModC01", as.integer(n), as.integer(p),
             as.integer(nrow(phi)),
            as.double(rho),
            phi=as.double(phi))
  res$phi <- matrix(res$phi, ncol=3)
  # edges are ranked according to the third empirical correlation 0-1
  res$phi <- res$phi[order(res$phi[,3],decreasing = TRUE),]
  dimnames(res$phi)[[2]] <- c("a","Vois(a)","lambda")
  return(res$phi)
}
