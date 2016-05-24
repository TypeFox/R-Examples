#############################################################
#
#	localdepth1Dsimplicialsimilarity function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: September, 2, 2008
#	Version: 0.1
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth1Dsimplicialsimilarity <- function(x, y, tau, use) {
  x <- as.vector(x)
  y <- as.vector(y)
  nx <- length(x)
  ny <- length(y)

  if (use=='diameter') nuse <- 1
  if (use=='volume') nuse <- 2
  if (use=='spherical') nuse <- 3
  
  ## number of couples
  nc <- choose(nx, 2)
  result <- .Fortran("LD1DSS",
    as.double(x),
    as.double(y),
    as.integer(nx), 
    as.integer(ny),
    as.integer(nc),
    as.double(tau),
    as.integer(nuse),               
    localdepth=mat.or.vec(ny,ny),               
    depth=mat.or.vec(ny,ny),
    diameters=double(nc),
    PACKAGE = "localdepth")
  result[[1]] <- result[[2]] <- result[[3]] <- result[[4]] <- NULL
  result[[1]] <- result[[2]] <- result[[3]] <- NULL
  result$areas <- result$diameters
  result$call <- match.call()
  result$num <- nc
  class(result) <- 'localdepth.similarity'
  return(result)
}
