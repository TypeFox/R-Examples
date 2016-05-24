#############################################################
#
#	quantile.linear2D function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: August, 26, 2013
#	Version: 0.1-1
#
#	Copyright (C) 2013 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

quantile.linear2D <- function(x, probs, nsamp='all', all=FALSE) {
  if (is.vector(x))
    x <- matrix(x, ncol=1)  
  x <- as.matrix(x)
  nrx <- nrow(x)
  if (nrx < 3) stop('x must have at least', 3, 'rows')
  if (is.character(nsamp)) {
    if (nsamp=='all') {
      nt <- choose(nrx, 3)
      nsamp <- FALSE
    } else {
      stop("if 'nsamp' is character then it must be equal to 'all'")
    }
  } else {
    if (nsamp < 1) {
      stop("'nsamp' must be greater than 0")
    } else { 
      nt <- nsamp
      nsamp <- TRUE
    }
  }
  if (nt > .Machine$integer.max/3)
    nt <- .Machine$integer.max/3

  if (nsamp) {
    result <- .Fortran("lldmcd2D",
      as.matrix(x),
      as.integer(nrx),
      as.integer(nt),
      result = double(nt*3),
      PACKAGE = "localdepth")$result
  } else {
    result <- .Fortran("lldmcd2D",
      as.matrix(x),
      as.integer(nrx),
      as.integer(nt),
      result = double(nt*3),
      PACKAGE = "localdepth")$result
  }
  res <- quantile(result, probs)
  if (all) {
     res <- list(quantile=res, stats=result, call=match.call())
  }
  class(res) <- 'quantile.localdepth'  
  return(res)
}

