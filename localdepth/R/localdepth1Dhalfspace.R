#############################################################
#
#	localdepth1Dhalfspace function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: March, 25, 2009
#	Version: 0.1
#
#	Copyright (C) 2009 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth1Dhalfspace <- function(x, y=NULL, tau) {
  x <- as.vector(x)
  y <- as.vector(y)
  nx <- length(x)
  ny <- length(y)
  result <- .Fortran("LD1DT",
    as.double(x),
    as.double(y),
    as.integer(nx), 
    as.integer(ny),
    as.double(tau),
    localdepth=double(ny),                   
    depth=double(ny),
    PACKAGE = "localdepth")
  result[[1]] <- result[[2]] <- result[[3]] <- result[[4]] <- result[[5]] <- NULL
  return(result)
}
