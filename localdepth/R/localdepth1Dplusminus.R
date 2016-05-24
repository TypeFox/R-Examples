#############################################################
#
#	localdepth1Dplusminus function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: November, 16, 2011
#	Version: 0.1
#
#	Copyright (C) 2011 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth1Dplusminus <- function(x, y=NULL, tau) {
  x <- as.vector(x)
  y <- as.vector(y)
  nx <- length(x)
  ny <- length(y)
  result <- .Fortran("LD1PM",
    as.double(x),
    as.double(y),
    as.integer(nx), 
    as.integer(ny),
    as.double(tau),
    poslocaldepth=double(ny),
    neglocaldepth=double(ny),
    posdepth=double(ny),
    negdepth=double(ny),
    PACKAGE = "localdepth")
  result[[1]] <- result[[2]] <- result[[3]] <- result[[4]] <- result[[5]] <- NULL
  return(result)
}
