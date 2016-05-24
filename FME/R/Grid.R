## -----------------------------------------------------------------------------
## Generates Parameters arranged on a Grid
## -----------------------------------------------------------------------------

Grid <- function(parRange, num) {

  npar  <- nrow(parRange)
  ngrid <- trunc(num^(1/npar))  # recalculate number of grid cells per parameter
  pargrid <- list()
  for (i in 1:npar) {
    pr   <- unlist(parRange[i,])
    pval <- seq(from = pr[1], to = pr[2], length = ngrid)
    pargrid[[i]] <- pval
  }
  pg <- as.matrix(expand.grid(pargrid))
  colnames(pg) <- rownames(parRange)
  return(pg)
}
