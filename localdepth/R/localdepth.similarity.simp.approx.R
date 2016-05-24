#############################################################
#
#	localdepth.similarity.simp.approx function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: November, 07, 2008
#	Version: 0.1-3
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth.similarity.simp.approx <- function(x, y=NULL, tau, use=c('volume', 'diameter'), nsamp='all', nmax=1, tol=0) {
  if (is.null(y))
    y <- x
  if (is.vector(x))
    x <- matrix(x, nrow=1)
  if (is.vector(y))
    y <- matrix(y, nrow=1)

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
  if (use=='diameter') use=0
  if (use=='volume') use=1
  
  if (is.character(nsamp) && nsamp=='all') {
    z <- .Fortran("ldssa",
      as.matrix(x), 
      as.matrix(y),
      as.double(tau),
      as.integer(nc),
      as.integer(nt),
      as.integer(nx),
      as.integer(ny),
      as.integer(use),
      as.double(tol),                  
      depth=mat.or.vec(ny,ny),
      localdepth=mat.or.vec(ny,ny),
      napprox=integer(ny),
      dapprox=double(1),
      PACKAGE = "localdepth")
  } else if (is.numeric(nsamp)){
    nt <- round(nt*nmax)
    z <- .Fortran("ldssaa",
      as.matrix(x), 
      as.matrix(y),
      as.double(tau),
      as.integer(nc),
      as.integer(nt),
      as.integer(nsamp),
      as.integer(nx),
      as.integer(ny),
      as.integer(use),
      as.double(tol),
      depth=mat.or.vec(ny,ny),
      localdepth=mat.or.vec(ny,ny),                  
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
  class(result) <- 'localdepth.similarity'
  return(result)
}
