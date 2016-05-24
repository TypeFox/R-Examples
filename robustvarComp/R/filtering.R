univariate.filter <- function(x, alpha=0.99) {
  pfnorm <- function(x) {
    2*pnorm(x) - 1
  }
  nx <- nrow(x)
  nc <- ncol(x)
  w <- matrix(FALSE, nx, nc)
  cc <- qnorm((1+alpha)/2)
  F <- ((1:nx)-1)/nx
  mx <- apply(x, 2, median)
  sx <- apply(x, 2, mad)
  z <- abs(scale(x, center=mx, scale=sx))
  zz <- apply(z, 2, sort)
  i0 <- apply(zz, 2, function(y) sum(y < cc))
  for (j in 1:nc) {
    ez <- pfnorm(zz[,j])
    d <- max(max(0,ez[(i0[j]+1):nx]-F[(i0[j]+1):nx]))
    if (length(zzz <- zz[ez >= 1 - d,j])) {
      tt <- min(zzz)      
      w[,j] <- z[,j] >= tt
    }
  }
  x[w] <- NA
  return(x)
}
