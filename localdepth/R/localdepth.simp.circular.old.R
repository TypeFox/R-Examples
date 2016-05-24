#############################################################
#
#	localdepth.simp.circular function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: October, 19, 2011
#	Version: 0.2
#
#	Copyright (C) 2011 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth.simp.circular.old <- function(x, y=NULL, tau, use=c('diameters', 'areas', 'spherical')) {
  use <- match.arg(use)
  if (is.null(y))
    y <- x
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  y <- conversion.circular(y, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(y, "class") <- attr(y, "circularp") <-  NULL
  if (!is.vector(x)) stop('y must be a vector')
  if (!is.vector(y)) stop('y must be a vector')
  nx <- length(x)
  if (nx < 2) stop('x must have at least length 2')
  ny <- length(y)
  
  ## number of couples
  nc <- choose(nx, 2)
  diameter <- rep(0, nc)
  depth <- ldepth <- rep(0, ny)
  ind1 <- 0
  for (i1 in 1:(nx-1)) {
    for (i2 in (i1+1):nx) {
      x1 <- x[i1]
      x2 <- x[i2]
      ind1<-ind1+1
      d21 <- abs((x2-x1)%%(2*pi))
      d12 <- abs((x1-x2)%%(2*pi))
      amb <- FALSE
      if (d21 == d12 & d12 == pi) {
        diameter[ind1] <- d21
        amb <- TRUE
      } else if (d21 < d12) {
        diameter[ind1] <- d21
      } else {
        diameter[ind1] <- d12
        z <- x1
        x1 <- x2
        x2 <- z
      }
      for (ind2 in 1:ny) {
        if (amb) {
          depth[ind2] <- depth[ind2]+1/2
          if (tau >= pi) ldepth[ind2] <- ldepth[ind2]+1/2
        } else { 
          if (abs((y[ind2]-x1)%%(2*pi)) <= diameter[ind1]) {
             depth[ind2] <- depth[ind2]+1
             if (use=='diamters' || use=='areas') {
               if (diameter[ind1] <= tau) {
                 ldepth[ind2] <- ldepth[ind2]+1
               }
             } else {
               spherical <- max(abs((y[ind2]-x1)%%(2*pi)), abs((x2-y[ind2])%%(2*pi)))
               if (spherical <= tau) {
                 ldepth[ind2] <- ldepth[ind2]+1
               }
             }
          }
        }
      }
    }
  }
  
  result <- list()
  result$num <- nc
  result$localdepth <- ldepth/nc
  result$depth <- depth/nc
  result$call <- match.call()
  result$tau <- tau
  result$use <- use
  result$diameters
  result$areas
  result$x <- x
  result$y <- y
  class(result) <- 'localdepth.simp.circular'
  return(result)
}
