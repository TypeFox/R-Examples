#############################################################
#
#	localdepth.halfspace function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: March, 25, 2009
#	Version: 0.1
#
#	Copyright (C) 2009 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth.halfspace <- function(x, y=NULL, tau, use=c('volume', 'diameter')) {
  use <- match.arg(use)
  if (is.null(y))
    y <- x
  if (is.vector(x)) {
    nx <- length(x)
    if (nx < 2) stop('x must have at least length 2')
    if (!is.vector(y)) stop('y must be a vector')
    result <- localdepth1Dhalfspace(x,y,tau)
  } else {
      stop('Not implemented for dimension different of 1')
  }

  result$localdepth <- result$localdepth/nx
  result$depth <- result$depth/nx
  result$max.localdepth <- max(result$localdepth)
  result$max.depth <- max(result$depth)  
  result$num <- c(nx, nx)
  result$call <- match.call()
  result$tau <- tau
  result$use <- use
  result$x <- x
  result$y <- y
  result$type <- 'exact'
  result$nsamp <- 'all'
  result$method <- 'halfspace'
  class(result) <- 'localdepth'
  return(result)
}

