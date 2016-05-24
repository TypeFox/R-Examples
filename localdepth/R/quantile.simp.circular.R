#############################################################
#
#	quantile.simp.circular function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: August, 26, 2013
#	Version: 0.3-1
#
#	Copyright (C) 2013 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

quantile.simp.circular <- function(x, probs, all=FALSE, ...) {
  nx <- length(x)
  diameters <- rep(0, choose(nx, 2))
  if (is.circular(x)) {
    xcp <- circularp(x)
  } else {
    xcp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
  }  
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  if (nx < 2) stop('x must have at least length 2')
  diameters <- .C("circdiam", x = as.double(x), nx = as.integer(nx),
               diameters = as.double(diameters),
               DUP = FALSE, NAOK = FALSE, PACKAGE = "localdepth")$diameters
  res <- quantile(diameters, probs, ...)
  res <- conversion.circular(circular(res), xcp$units, xcp$type, xcp$template, xcp$modulo, xcp$zero, xcp$rotation)
  if (all) {
    diameters <- conversion.circular(circular(diameters), xcp$units, xcp$type, xcp$template, xcp$modulo, xcp$zero, xcp$rotation)
     res <- list(quantile=res, stats=diameters, call=match.call())
  }
  class(res) <- 'quantile.localdepth'
  return(res)
}

