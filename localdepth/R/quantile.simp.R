#############################################################
#
#	quantile.simp function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: August, 26, 2013
#	Version: 0.5-2
#
#	Copyright (C) 2013 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

quantile.simp <- function(x, probs, use=c('volume', 'diameter'), nsamp='all', all=FALSE, ...) {
  use <- match.arg(use)
##  norm <- function(x) sqrt(t(x)%*%x)
##  area <- function(x) sqrt(sum(x)*(sum(x)/2-x[1])*(sum(x)/2-x[2])*(sum(x)/2-x[3])/2)
  x <- as.matrix(x)
  nx <- nrow(x)
  nc <- ncol(x)
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
          result <- .Fortran("lddiamsa",
                    as.double(x),
                    as.integer(nt),
                    as.integer(nc),
                    as.integer(nx),
                    result = double(nt),
                   PACKAGE = "localdepth")$result
        } else {
          result <- .Fortran("ldareasa",
                    as.double(x),
                    as.integer(nt),
                    as.integer(nc),
                    as.integer(nx),
                    result = double(nt),
                   PACKAGE = "localdepth")$result
        }
  } else {
    if (is.matrix(x) | is.data.frame(x)) {
      if (ncol(x)==2) {
        result <- rep(0, nt)
        if (use=='diameter') {
          result <- .C("twoDdiam", x = as.double(x[,1]), y = as.double(x[,2]),
               nx = as.integer(nx), result = as.double(result),
               DUP = FALSE, NAOK = FALSE, PACKAGE = "localdepth")$result
        } else {
          result <- .C("twoDarea", x = as.double(x[,1]), y = as.double(x[,2]),
               nx = as.integer(nx), result = as.double(result),
               DUP = FALSE, NAOK = FALSE, PACKAGE = "localdepth")$result
        }
      } else {
        if (use=='diameter') {
          result <- .Fortran("lddiams",
                    as.double(x),
                    as.integer(nt),
                    as.integer(nc),
                    as.integer(nx),
                    result = double(nt),
                   PACKAGE = "localdepth")$result
        } else {
          result <- .Fortran("ldareas",
                    as.double(x),
                    as.integer(nt),
                    as.integer(nc),
                    as.integer(nx),
                    result = double(nt),
                   PACKAGE = "localdepth")$result
        }
      }
    } else if (is.vector(x)) {
      nx <- length(x)
      if (nx < 2) stop('x must have at least length 2')
      nc <- choose(nx, 2)
      result <- rep(0, nc)
      result <- .C("oneDdiam", x = as.double(x), nx = as.integer(nx),
               result = as.double(result),
               DUP = FALSE, NAOK = FALSE, PACKAGE = "localdepth")$result
    }
  }
  res <- quantile(result, probs, ...)
  if (all) {
     res <- list(quantile=res, stats=result, call=match.call())
  }
  class(res) <- 'quantile.localdepth'  
  return(res)
}

