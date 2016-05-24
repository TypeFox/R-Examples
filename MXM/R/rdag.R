rdag <- function(n, p, s, a = 0, m, A = NULL, seed = FALSE) {
  ## n in the sample size
  ## p is the number of (nodes or) variables
  ## s is the success probability of the binomial distribution
  ## a is the percentage of outliers, is set to zero by default
  ## a number between 0 and 1
  ## m is the mean vector which is used only if you want outliers, i.e. if a > 0
  ## A is an adjancey matrix given by the user 
  
  if ( is.null(A) ) { ## no adjacency matrix is given
    if ( s > 1 || s < 0 )  s <- 0.5
    if ( a > 1 || a < 0 )  a <- 0
    if ( seed == TRUE )  set.seed(1234567)
    A <- matrix( numeric(p^2), nrow = p )
    nu <- 0.5 * p * (p - 1)
    A[ lower.tri(A) ] <- rbinom(nu, 1, s)
    A[ A == 1 ] <- runif( sum(A), 0.1, 1 )
  } else {
    A <- A
    p <- ncol(A)
  }
   
  Ip <- diag(p)
  sigma <- solve( Ip - A )
  sigma <- crossprod( sigma ) 
  nout <- 0
  if ( seed == TRUE )  set.seed(1234567)
  
  if (a > 0) {
    y <- MASS::mvrnorm( n - nout, numeric(p), sigma)
    nout <- round( a * n )
    yout <- MASS::mvrnorm(nout, m, sigma)
    x <- rbind(y, yout)  
  } else {
    x <- MASS::mvrnorm(n, numeric(p), sigma)
  }
  
  B <- A
  B[ B > 0 ] <- 2
  ind <- which( t(B) == 2 )
  B[ind] <- 3
  
  colnames(x) <- paste("X", 1:p, sep = "")
  list(nout = nout, G = t( B ), A = A, x = x)
}