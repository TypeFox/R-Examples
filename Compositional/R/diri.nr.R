################################
#### Dirichlet distribution parameters 
#### via Newton-Raphson
#### Tsagris Michail 8/2015  
#### mtsagris@yahoo.gr
#### References: Estimating a Dirichlet distribution (2012)
#### Thomas P. Minka 
#### http://research.microsoft.com/en-us/um/people/minka/papers/dirichlet/minka-dirichlet.pdf
################################

diri.nr <- function(x) {
  ## x is compositional data

  runtime <- proc.time()
  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x/rowSums(x)  ## makes sure x is compositional data
  n <- nrow(x)  ## the sample size 
  p <- ncol(x)  ## dimensionality
  ell <- NULL
  m <- colMeans(x)
  zx <- log(x)
  
  lm <- log(m)
  down <-  - sum(  m * colMeans( log( x ) - rep(lm, rep(n, p ) ) ) )
  s <- 0.5 * (p - 1) / down  ## initial value for precision
  a <- s * m  ## initial values
  ell[1] <- n * lgamma( sum(a) ) - n * sum( lgamma(a) ) +
  sum( zx %*% (a- 1) ) 
  gm <- colSums(zx)
  z <- n * digamma( sum(a) )
  g <- z - n * digamma(a) + gm
  qk <-  - n * trigamma(a)
  b <- ( sum(g / qk) ) / ( 1/z - sum(1 / qk) )
  a <- a - (g - b)/qk
  ell[2] <- n * lgamma( sum(a) ) - n * sum( lgamma(a) ) +
  sum( zx %*% (a - 1) ) 
  i <- 2

  while( abs(ell[i] - ell[i - 1]) > 1e-07 ) {
    i <- i + 1
    z <- n * digamma( sum(a) )
    g <- z - n * digamma(a) + gm
    qk <-  - n * trigamma(a)
    b <- ( sum(g / qk) ) / ( 1/z - sum(1 / qk) )
    a <- a - (g - b) / qk
    ell[i] <- n * lgamma( sum(a) ) - n * sum( lgamma(a) ) +
    sum( zx %*% (a - 1) ) 
  }

  if ( is.null(colnames(x)) ) {
    names(a) <- paste("X", 1:p, sep = "")
  } else  names(a) <- colnames(x)
  
  runtime <- proc.time() - runtime
  list(iter = i, loglik = ell[i], param = a, runtime = runtime)
}