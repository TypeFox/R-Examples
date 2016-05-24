#############################################################
#                                                           #
#   ticks.circular function                                 #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: November, 13, 2008                                #
#   Version: 0.4-1                                          #
#                                                           #
#   Copyright (C) 2008 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
ticks.circular <- function(x, template=c("none", "geographics"), zero=NULL, rotation=NULL, tcl=0.025, col=NULL, ...) {
  template <- match.arg(template)
  if (is.null(col)) col <- par("col")
  xcircularp <- attr(as.circular(x), "circularp")
  type <- xcircularp$type
  if (type=='directions')
    x <- 2*x
  if (template=="geographics") {
    zero <- pi/2
    rotation <- "clock"
  } else {
    if (is.null(zero))
      zero <- xcircularp$zero
    if (is.null(rotation))
      rotation <- xcircularp$rotation
  }  
  x <- conversion.circular(x, units="radians")  
  attr(x, "circularp") <- attr(x, "class") <- NULL
  if (rotation=="clock")
    x <- -x
  x <- x+zero
  x <- x%%(2*pi)
  TicksCircularRad(x, tcl, col, ...)
}

TicksCircularRad <- function(x, tcl, col, ...) {
    r <- 1+tcl*c(-1/2,1/2)
    z <- cos(x)
    y <- sin(x)
    for (i in 1:length(x)) {
       lines.default(z[i]*r, y[i]*r, col=col, ...)
    }
}

