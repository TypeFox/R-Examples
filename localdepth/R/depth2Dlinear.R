#############################################################
#
#	depth2Dlinear function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: November, 20, 2008
#	Version: 0.1
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

depth2Dlinear <- function(x, y=NULL, nsamp='all', nmax=1) {
  if (is.null(y))
    y <- x
  if (is.vector(x))
    x <- matrix(x, ncol=1)
  if (is.vector(y))
    y <- matrix(y, ncol=1)
  x <- as.matrix(x)
  y <- as.matrix(y)
  nrx <- nrow(x)
  if (nrx < 3) stop('x must have at least', 3, 'rows')  
  nry <- nrow(y)
  result <- list()
  depth <- rep(0,nry)
  nt <- choose(nrx, 3)
  if (nt > .Machine$integer.max/3)
    nt <- .Machine$integer.max/3
  nt <- nt*nmax
  if (is.numeric(nsamp) && nsamp <= 0) stop("the argument 'nsamp' must be positive")
  if (is.numeric(nsamp) && nsamp > nt) {
      warning("Since 'nsamp' is greater than the number of simplex the 'all' method is used")
      nsamp <- 'all'
  }

  if (is.character(nsamp) && nsamp=='all') {  
    z <- .Fortran("dlal2D",
      as.matrix(x),
      as.matrix(y),
      as.integer(nrx),
      as.integer(nry),
      as.integer(nsamp),
      depth = double(nry),
      PACKAGE = "localdepth")
  } else {
    z <- .Fortran("dlmc2D",
      as.matrix(x),
      as.matrix(y),
      as.integer(nrx),
      as.integer(nry),
      as.integer(nsamp),
      depth = double(nry),
      PACKAGE = "localdepth")
  }

  result$depth <- z$depth/nt
  result$max.depth <- max(result$depth)
  result$num <- c(nt,nt)    
  result$call <- match.call()
  result$x <- x
  result$y <- y
  result$type <- 'exact'
  result$nsamp <- nsamp
  result$method <- 'linear'
  class(result) <- 'localdepth'
  return(result)
}

