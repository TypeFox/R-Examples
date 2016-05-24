#############################################################
#
#	depth2Dsigma function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: July, 05, 2013
#	Version: 0.1
#
#	Copyright (C) 2013 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

depth2Dsigma <- function(x, y, sigma) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  nx <- nrow(x)
  ny <- nrow(y)
  result <- list()
  ## numero di terne
  nt <- choose(nx, 3)
  depth <- rep(0,ny)
  res <- .C("d2Dsigma",
              x1 = as.double(x[,1]),
              x2 = as.double(x[,2]),
              y1 = as.double(y[,1]),
              y2 = as.double(y[,2]),
              sigma =as.double(sigma),
              nx = as.integer(nx),
              ny = as.integer(ny),
              depth = as.double(depth),
              DUP = FALSE, NAOK = FALSE, PACKAGE = "localdepth")
  return(res$depth/nt)
}

