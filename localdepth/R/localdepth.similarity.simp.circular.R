#############################################################
#
#	localdepth.similarity.simp.circular function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: August, 28, 2008
#	Version: 0.1
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth.similarity.simp.circular <- function(x, y=NULL, tau, use=c('volume', 'diameter', 'spherical')) {
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

  if (use=='diameter') nuse <- 1
  if (use=='volume') nuse <- 2
  if (use=='spherical') nuse <- 3
  
  ## number of couples
  nc <- choose(nx, 2)
  diameters <- rep(0, nc)
  depth <- rep(0, ny*ny)
  localdepth <- rep(0, ny*ny)
  
  res <- .C("ldcircsimpsim", x = as.double(x), y = as.double(y),
            nx = as.integer(nx), ny = as.integer(ny), tau = as.double(tau),
            nuse = as.integer(nuse), depth = as.double(depth),
            localdepth = as.double(localdepth), diameters = as.double(diameters),
            DUP = TRUE, NAOK = FALSE, PACKAGE = "localdepth")

  depth <- matrix(res$depth, ny, ny, byrow=FALSE)
  depth[lower.tri(depth)] <- t(depth)[lower.tri(depth)]
  localdepth <- matrix(res$localdepth, ny, ny, byrow=FALSE)
  localdepth[lower.tri(localdepth)] <- t(localdepth)[lower.tri(localdepth)]

  result <- list()
  result$localdepth <- localdepth/nc
  result$depth <- depth/nc
  result$max.localdepth <- max(result$localdepth)
  result$max.depth <- max(result$depth)
  result$num <- nc
  result$call <- match.call()
  result$tau <- tau
  result$use <- use
  result$type <- 'exact'
  result$nsamp <- 'all'
  result$method <- 'simplicial'  
  class(result) <- 'localdepth.similarity'
  return(result)
}
