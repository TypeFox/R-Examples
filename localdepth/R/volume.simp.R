#############################################################
#
#	volume.simp function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: June, 23, 2008
#	Version: 0.1
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

volume.simp <- function(x) {
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  if (nr!=(nc+1)) stop("the matrix 'x' must be (n+1)*n dimension")
  x <- abs(det(cbind(rep(1, nr), x)))/factorial(nc)
  return(x)
}
