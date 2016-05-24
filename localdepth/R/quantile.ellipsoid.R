#############################################################
#
#	quantile.ellipsoid function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: August, 26, 2013
#	Version: 0.1-4
#
#	Copyright (C) 2013 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

quantile.ellipsoid <- function(x, probs, use=c('volume', 'diameter'), nsamp='all', all=FALSE, dimension=NULL, ...) {
  use <- match.arg(use)
  if (is.vector(x))
    x <- matrix(x, ncol=1)  
  x <- as.matrix(x)
  nx <- nrow(x)
  nc <- ncol(x)
  if (is.null(dimension))
    dimension <- nc  
  if (nx < nc+1) stop('x must have at least', nc+1, 'rows')
  if (is.character(nsamp)) {
    if (nsamp=='all') {
      nt <- choose(nx, nc+1)
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
          result <- .Fortran("eldiamsa",
                    as.double(x),
                    as.integer(nt),
                    as.integer(nc),
                    as.integer(nx),
                    as.double(dimension),
                    result = double(nt),
                   PACKAGE = "localdepth")$result
        } else {
          result <- .Fortran("elareasa",
                    as.double(x),
                    as.integer(nt),
                    as.integer(nc),
                    as.integer(nx),
                    as.double(dimension),
                    result = double(nt),
                   PACKAGE = "localdepth")$result
        }
  } else {
        if (use=='diameter') {
          result <- .Fortran("eldiams",
                    as.double(x),
                    as.integer(nt),
                    as.integer(nc),
                    as.integer(nx),
                    as.double(dimension),
                    result = double(nt),
                   PACKAGE = "localdepth")$result
        } else {
          result <- .Fortran("elareas",
                    as.double(x),
                    as.integer(nt),
                    as.integer(nc),
                    as.integer(nx),
                    as.double(dimension),
                    result = double(nt),
                   PACKAGE = "localdepth")$result
        }
  }
  res <- quantile(result, probs, ...)
  if (all) {
     res <- list(quantile=res, stats=result, call=match.call())
  }
  class(res) <- 'quantile.localdepth'
  return(res)
}

