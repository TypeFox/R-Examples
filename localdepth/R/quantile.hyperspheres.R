#############################################################
#
#	quantile.hyperspheres function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: August, 26, 2013
#	Version: 0.1-1
#
#	Copyright (C) 2013 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

quantile.hyperspheres <- function(x, probs, use=c('volume', 'diameter'), nsamp='all', size=FALSE, ...) {
  use <- match.arg(use)
  if (is.vector(x))
    x <- matrix(x, ncol=1)  
  x <- as.matrix(x)
  nx <- nrow(x)
  nc <- ncol(x)
  if (nx < nc) stop('x must have at least', nc, 'rows')
  if (use=='volume' & nc!=3)
    stop("The option use='volume' is available only on the sphere")  
  if (is.character(nsamp)) {
    if (nsamp=='all') {
      nt <- choose(nx, nc)
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
  
  if (nsamp) {
        if (use=='diameter') {
          result <- .Fortran("diamshsa",
                    as.double(x),
                    as.integer(nt),
                    as.integer(nc),
                    as.integer(nx),
                    result = double(nt),
                   PACKAGE = "localdepth")$result
        } else {
          result <- .Fortran("areashsa",
                    as.double(x),
                    as.integer(nt),
                    as.integer(nc),
                    as.integer(nx),
                    result = double(nt),
                   PACKAGE = "localdepth")$result
        }
  } else {
        if (use=='diameter') {
          result <- .Fortran("lddiamshs",
                    as.double(x),
                    as.integer(nt),
                    as.integer(nc),
                    as.integer(nx),
                    result = double(nt),
                   PACKAGE = "localdepth")$result
        } else {
          result <- .Fortran("ldareashs",
                    as.double(x),
                    as.integer(nt),
                    as.integer(nc),
                    as.integer(nx),
                    result = double(nt),
                   PACKAGE = "localdepth")$result
        }
  }
  res <- quantile(result, probs, ...)
  if (size) {
     res <- list(quantile=res, stats=result, call=match.call())
  }
  class(res) <- 'quantile.localdepth'
  return(res)
}
