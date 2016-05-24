#############################################################
#
#	localdepth1Dfprime function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: November, 16, 2011
#	Version: 0.1
#
#	Copyright (C) 2011 Claudio Agostinelli, Juan Francisco Rosco Nieves and Mario Romanazzi
#
#############################################################

localdepth1Dfprime <- function(x, y=NULL, tau) {
  x <- as.vector(x)
  y <- as.vector(y)
  nx <- length(x)
  ny <- length(y)
  result <- .Fortran("LD1FP",
    as.double(x),
    as.double(y),
    as.integer(nx), 
    as.integer(ny),
    as.double(tau),
    fprime=double(ny),
    PACKAGE = "localdepth")
  result[[1]] <- result[[2]] <- result[[3]] <- result[[4]] <- result[[5]] <- NULL
  return(result$fprime/ny)
}
