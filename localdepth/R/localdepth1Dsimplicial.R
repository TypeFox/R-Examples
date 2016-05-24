#############################################################
#
#	localdepth1Dsimplicial function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: February, 10, 2011
#	Version: 0.2-1
#
#	Copyright (C) 2011 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth1Dsimplicial <- function(x, y, tau, use) {
  x <- as.vector(x)
  y <- as.vector(y)
  nx <- length(x)
  ny <- length(y)

  if (use=='diameter') nuse <- 1
  if (use=='volume') nuse <- 2
  if (use=='spherical') nuse <- 3
  
  ## number of couples
  nc <- choose(nx, 2)
  result <- .Fortran("LD1DS",
    as.double(x),
    as.double(y),
    as.integer(nx), 
    as.integer(ny),
    as.integer(nc),
    as.double(tau),
    as.integer(nuse),               
    localdepth=double(ny),               
    depth=double(ny),
    diameters=double(nc),
    PACKAGE = "localdepth")
  result[[1]] <- result[[2]] <- result[[3]] <- result[[4]] <- NULL
  result[[1]] <- result[[2]] <- result[[3]] <- NULL
  result$diameters <- NULL
  ##result$areas <- result$diameters
  result$num <- nc
  return(result)
}

