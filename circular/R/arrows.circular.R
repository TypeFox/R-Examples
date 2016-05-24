#############################################################
#                                                           #
#   arrows.circular function                                #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: August, 14, 2010                                  #
#   Version: 0.3                                            #
#                                                           #
#   Copyright (C) 2010 Claudio Agostinelli                  #
#                                                           #
#############################################################
# patche suggests by Peter Cowan (pdc)
# [#193] fix for plotting many arrows with one call to arrows.circular()

arrows.circular <- function(x, y=NULL, x0=0, y0=0, na.rm=FALSE, shrink=1, plot.info=NULL, zero=NULL, rotation=NULL, ...) {
  if (na.rm)
    x <- x[!is.na(x)]
  if (length(x)==0) {
    warning("No observations (at least after removing missing values)")
    return(NULL)
  }
  xcircularp <- attr(as.circular(x), "circularp")
  if (is.null(plot.info)) {
    if (is.null(zero))
      zero <- xcircularp$zero
    if (is.null(rotation))
      rotation <- xcircularp$rotation
  } else {
    zero <- plot.info$zero
    rotation <- plot.info$rotation
  }
  x <- conversion.circular(x, units="radians")
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  if (rotation=="clock")
    x <- -x
  x <- x+zero
  x <- x%%(2*pi)
  x <- as.vector(x)
  if (is.null(y))
    y <- rep(1, length(x))
  y <- as.vector(y)
  if (length(y)!=length(x))
    stop("'y' must have the same length of 'x'")
  y <- y*shrink
  if (length(x0)!=length(x))
    x0 <- rep(x0, length(x))
  if (length(y0)!=length(x))
    y0 <- rep(y0, length(x))
  x1 <- x0 + y*cos(x)
  y1 <- y0 + y*sin(x)
  arrows(x0, y0, x1, y1, ...)  
}
