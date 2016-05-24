################################
#### Symmetric Dirichlet distribution
#### Tsagris Michail 11/2013
#### mtsagris@yahoo.gr
################################

sym.test <- function(x) {
  ## x contains the data
  n <- nrow(x)  ## the sample size
  D <- ncol(x)  ## the dimensionality of the data
  zx <- log(x)
   sym <- function(a, zx) {
    n * lgamma(D * a) - n * D * lgamma(a) +
    sum( zx * (a - 1) )
   }
  t0 <- optimize(sym, c(0, 1000), zx = zx, maximum = TRUE)
  t1 <- diri.nr(x)
  a0 <- t0$maximum
  a1 <- t1$param
  h1 <- t1$loglik
  h0 <- as.numeric(t0$objective)
  test <- 2 * (h1 - h0)
  pvalue <- 1 - pchisq(test, D - 1)

  if ( is.null(colnames(x)) ) {
    names(a1) <- paste("X", 1:D, sep = "")
  } else  names(a1) <- colnames(x)
  res <- c(h1, h0, test, pvalue, D - 1)
  names(res) <- c('loglik1', 'loglik0', 'test', 'pvalue', 'df')
  list(est.par = a1, one.par = a0, res = res )
}
