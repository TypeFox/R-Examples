#############################################################
#
#	localdepth.simp.circular function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: July, 19, 2011
#	Version: 0.4
#
#	Copyright (C) 2011 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth.simp.circular <- function(x, y=NULL, tau, use=c('volume', 'diameter', 'spherical')) {
  use <- match.arg(use)
  if (is.null(y))
    y <- x
  if (is.circular(x)) {
    xcp <- circularp(x)
  } else {
    xcp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
  }
  if (is.circular(y)) {
    ycp <- circularp(y)
  } else {
    ycp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
  }
  if (is.circular(tau)) {
    tcp <- circularp(tau)
  } else {
    tcp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
  }
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  y <- conversion.circular(y, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(y, "class") <- attr(y, "circularp") <-  NULL
  tau <- conversion.circular(tau, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(tau, "class") <- attr(tau, "circularp") <-  NULL
  if (!is.vector(x)) stop('y must be a vector')
  if (!is.vector(y)) stop('y must be a vector')
  nx <- length(x)
  if (nx < 2) stop('x must have at least length 2')
  ny <- length(y)
  if (ny < 1) stop('y must have at least length 1')
  if (length(tau)!=1) stop('tau must have length 1')

  if (use=='diameter') nuse <- 1
  if (use=='volume') nuse <- 2
  if (use=='spherical') nuse <- 3
  
  ## number of couples
  nc <- choose(nx, 2)
  diameters <- rep(0, nc)
  depth <- rep(0, ny)
  localdepth <- rep(0, ny)
  
  res <- .C("ldcircsimp", x = as.double(x), y = as.double(y),
            nx = as.integer(nx), ny = as.integer(ny), tau = as.double(tau),
            nuse = as.integer(nuse), depth = as.double(depth),
            localdepth = as.double(localdepth), diameters = as.double(diameters),
            DUP = FALSE, NAOK = FALSE, PACKAGE = "localdepth")
  
  result <- list()
  result$localdepth <- res$localdepth/nc
  result$depth <- res$depth/nc
  result$max.localdepth <- max(result$localdepth)
  result$max.depth <- max(result$depth)    
  result$num <- nc
  result$call <- match.call()
  result$tau <- conversion.circular(circular(tau), tcp$units, tcp$type, tcp$template, tcp$modulo, tcp$zero, tcp$rotation)
  result$use <- use
  result$x <- conversion.circular(circular(x), xcp$units, xcp$type, xcp$template, xcp$modulo, xcp$zero, xcp$rotation)
  result$y <- conversion.circular(circular(y), ycp$units, ycp$type, ycp$template, ycp$modulo, ycp$zero, ycp$rotation)
  result$type <- 'exact'
  result$nsamp <- 'all'
  result$method <- 'simplicial'  
  class(result) <- 'localdepth'
  return(result)
}
