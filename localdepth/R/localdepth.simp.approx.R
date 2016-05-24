#############################################################
#
#	localdepth.simp.approx function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: November, 7, 2008
#	Version: 0.1-6
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth.simp.approx <- function(x, y=NULL, tau, use=c('volume', 'diameter'), nsamp='all', nmax=1, tol=0) {
  if (is.null(y))
    y <- x
  if (is.vector(x))
    x <- matrix(x, ncol=1)
  if (is.vector(y))
    y <- matrix(y, ncol=1)
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  nx <- nrow(x)
  ny <- nrow(y)
  nc <- ncol(x)
  nt <- choose(nx, nc+1)
  if (nt > .Machine$integer.max)
    nt <- .Machine$integer.max

  use <- match.arg(use)

  if (is.numeric(nsamp) && nsamp <= 0) stop("the argument 'nsamp' must be positive")
  if (is.numeric(nsamp) && nsamp > nt) {
      warning("Since 'nsamp' is greater than the number of simplex the 'all' method is used")
      nsamp <- 'all'
  }
  if (use=='diameter') nuse=0
  if (use=='volume') nuse=1

#  if (est.time) {
#    est.time <- system.time(zest <- .Fortran("ldsaa",
#      as.matrix(x), 
#      as.matrix(y),
#      as.double(tau),
#      as.integer(nc),
#      as.integer(10000),
#      as.integer(10000),
#      as.integer(nx),
#      as.integer(ny),
#      as.integer(nuse),         
#      depth=double(ny),
#      localdepth=double(ny),                  
#      nd=double(1),
#      nld=double(1),
#      napprox=integer(ny),
#      dapprox=double(1),                  
#      PACKAGE = "localdepth"))
#    cat(est.time, '\n')
#    if (is.character(nsamp) && nsamp=='all') nsamp <- nt
#    cat('Estimated time: ', est.time[3]/10000*min(nt*nmax, nsamp)*5, '\n')
#  }
  
  if (is.character(nsamp) && nsamp=='all') {
    z <- .Fortran("ldsa",
      as.matrix(x), 
      as.matrix(y),
      as.double(tau),
      as.integer(nc),
      as.integer(nt),
      as.integer(nx),
      as.integer(ny),
      as.integer(nuse),
      as.double(tol),
      depth=double(ny),
      localdepth=double(ny),
      napprox=integer(ny),
      dapprox=double(1),
      PACKAGE = "localdepth")
  } else if (is.numeric(nsamp)){
    nt <- round(nt*nmax)
    z <- .Fortran("ldsaa",
      as.matrix(x), 
      as.matrix(y),
      as.double(tau),
      as.integer(nc),
      as.integer(nt),
      as.integer(nsamp),
      as.integer(nx),
      as.integer(ny),
      as.integer(nuse),
      as.double(tol),
      depth=double(ny),
      localdepth=double(ny),                  
      nd=double(1),
      nld=double(1),
      napprox=integer(ny),
      dapprox=double(1),                  
      PACKAGE = "localdepth") 
  } else {
    stop("the argument 'nsamp' must be either 'all' or a positive number")
  }
  result <- list()
  if (is.numeric(nsamp)) {
    result$localdepth <- z$localdepth/z$nd
    result$depth <- z$depth/z$nd
    result$max.localdepth <- max(result$localdepth)
    result$max.depth <- max(result$depth)  
    result$num <- c(z$nd,z$nld)
  } else {
    result$localdepth <- z$localdepth/nt
    result$depth <- z$depth/nt
    result$max.localdepth <- max(result$localdepth)
    result$max.depth <- max(result$depth)
    result$num <- c(nt,NA)
  }
  result$call <- match.call()
  result$tau <- tau
  result$use <- use
  result$tol <- tol  
  result$napprox <- z$napprox
  result$prob.approx <- z$dapprox
  result$x <- x
  result$y <- y
  result$type <- 'approx'
  result$nsamp <- nsamp
  result$method <- 'simplicial'
  class(result) <- 'localdepth'
  return(result)
}



### old implementation will be deleted very soon from the source code
localdepth.simp.approx.internal <- function(x, y) {
  p <- NCOL(x)
  depth <- rep(0, NROW(y))
  mx <- mean(as.data.frame(x))
  Prec <- solve(cov(x)*p/(p+1))
  y <- t(y)-mx
  mah <- function(x, Prec) {
    t(x)%*%Prec%*%x
  }
  mahy <- apply(X=y, MARGIN=2, FUN=mah, Prec=Prec)
  circum2 <- c(t(x[1,]-mx) %*% Prec %*% (x[1,]-mx))

### http://en.wikipedia.org/wiki/Simplex  
  volumesimplex <- sqrt(circum2)^(p+1)*sqrt(p+1)/factorial(p)*2^(p/2)
  volumecircum <- pi^(p/2)*((circum2)^(p/2))/gamma(p/2+1)
  volumeinner <- volumecircum/p^p

###  cat(volumeinner, volumesimplex, volumecircum, "\n")
  
  depth[mahy < circum2/(p^2)] <- 1
  depth[mahy >= circum2/(p^2) & mahy < circum2] <- (volumesimplex-volumeinner)/(volumecircum - volumeinner)
  return(depth)
}
