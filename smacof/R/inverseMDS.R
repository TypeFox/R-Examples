inverseMDS <- function (x) {
  x <- as.data.frame(x)
  n <- nrow(x)
  m <- ncol(x)
  x <- apply(x, 2, function (y) y - mean (y))
  nm <- n - (m + 1)
  kk <- cbind(1, x, matrix (rnorm (n * nm), n, nm))
  kperp <- as.matrix(qr.Q(qr(kk))[, -(1:(m + 1))])
  dd <- Euclid(x)
  k <- 1
  base <- matrix (0, n * (n - 1) / 2, nm * (nm + 1) / 2)
  for (i in 1 : nm) {
    for (j in 1 : i) {
      oo <- outer(kperp[, i], kperp[, j])
      if (j != i) {
        oo <- oo + t(oo)
      }
      base[, k] <- lower_triangle(dd * (1 - oo))
      k <- k + 1
    }
  }
  base <- cbind(lower_triangle(dd), base)
  baselist <- list()
  
  for (i in 1:ncol(base)) {
    dmat <- matrix(NA, n, n)
    dmat[lower.tri(dmat)] <- base[,1]
    baselist[[i]] <- as.dist(dmat)
    attr(baselist[[i]], "Labels") <- rownames(x)
  }  
  return (base = baselist)
}


