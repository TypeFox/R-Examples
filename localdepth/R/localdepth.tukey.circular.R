#############################################################
#
#	localdepth.tukey.circular function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: July, 19, 2011
#	Version: 0.2
#
#	Copyright (C) 2011 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

localdepth.tukey.circular <- function(x, y=NULL, tau, tol) {
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
  if (!is.vector(x)) stop('x must be a vector')
  if (!is.vector(y)) stop('y must be a vector')
  nx <- length(x)
  if (nx < 2) stop('x must have at least length 2')
  ny <- length(y)
  if (ny < 1) stop('y must have at least length 1')
  if (length(tau)!=1) stop('tau must have length 1')
    
  res <- .Fortran("ldtc",
    x = as.double(x),
    y = as.double(y),
    nx = as.integer(nx),
    ny = as.integer(ny),
    tau = as.double(tau),
    tol = as.double(tol),              
    depth = double(ny),
    localdepth = double(ny),
    PACKAGE = "localdepth")
  
  result <- list()
  result$localdepth <- res$localdepth/nx
  result$depth <- res$depth/nx
  result$max.localdepth <- max(result$localdepth)
  result$max.depth <- max(result$depth)    
  result$num <- nx
  result$call <- match.call()
  result$tau <- conversion.circular(circular(tau), tcp$units, tcp$type, tcp$template, tcp$modulo, tcp$zero, tcp$rotation)
  result$use <- NA
  result$x <- conversion.circular(circular(x), xcp$units, xcp$type, xcp$template, xcp$modulo, xcp$zero, xcp$rotation)
  result$y <- conversion.circular(circular(y), ycp$units, ycp$type, ycp$template, ycp$modulo, ycp$zero, ycp$rotation)
  result$type <- 'exact'
  result$nsamp <- 'all'
  result$method <- 'halfspace'  
  class(result) <- 'localdepth'
  return(result)
}
