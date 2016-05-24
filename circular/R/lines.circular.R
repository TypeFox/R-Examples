#############################################################
#                                                           #
#   lines.circular function                                 #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: March, 29, 2010                                   #
#   Version: 0.2-1                                          #
#                                                           #
#   Copyright (C) 2010 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
lines.circular <- function(x, y, join=FALSE, nosort=FALSE, offset=1, shrink=1, plot.info=NULL, zero=NULL, rotation=NULL, modulo=NULL, ...) {

   xcircularp <- attr(as.circular(x), "circularp")
#   type <- xcircularp$type
   if (is.null(modulo))
      modulo <- xcircularp$modulo
   if (is.null(plot.info)) {
      if (is.null(zero))
         zero <- xcircularp$zero
      if (is.null(rotation))
         rotation <- xcircularp$rotation
      next.points <- 0
   } else {
      zero <- plot.info$zero
      rotation <- plot.info$rotation
      next.points <- plot.info$next.points
   }

  
# Handling missing values
   ok <- complete.cases(x, y)
   x <- x[ok]
   y <- y[ok]

   if (length(x)) {
      x <- conversion.circular(x, units="radians", modulo=modulo)
      attr(x, "circularp") <- attr(x, "class") <- NULL
      attr(y, "circularp") <- attr(y, "class") <- NULL
      if (rotation=="clock")
         x <- -x
      x <- x+zero
###      x <- x%%(2*pi)
      ll <- LinesCircularRad(x, y, join, nosort, offset, shrink, ...) 
   }
   return(invisible(list(x=ll$x, y=ll$y, zero=zero, rotation=rotation, next.points=next.points)))
}

LinesCircularRad <- function(x, y, join=FALSE, nosort=FALSE, offset=1, shrink=1, ...) {
   n <- length(x)
   if (!nosort) {
      xorder <- order(x)
      x <- x[xorder]
      y <- y[xorder]
      spacings <- c(diff(x), x[1] - x[n] + 2*pi)
      pos <- which.max(spacings)[1]
      if (pos==n)
         xorder <- 1:n
      else
         xorder <- c((pos+1):n, 1:pos)
   } else {
         xorder <- 1:n
   
   }
   z <- (y/shrink+offset)*cos(x)
   w <- (y/shrink+offset)*sin(x)
   z <- z[xorder]
   w <- w[xorder]
   if (join) {
      z <- c(z, z[1])
      w <- c(w, w[1])
   }
   lines.default(x=z, y=w, ...)
   invisible(list(x=z, y=w))
}

