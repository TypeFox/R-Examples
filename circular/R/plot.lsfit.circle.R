
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   plot.lsfit.circle function                              #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: December, 15, 2004                                #
#   Version: 0.1                                            #
#                                                           #
#   Copyright (C) 2004 Claudio Agostinelli                  #
#                                                           #
#############################################################


plot.lsfit.circle <- function(x, add=FALSE, main=NULL, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, uin, tol=0.04, plus.cex=1, ...){

   xx <- x
   x <- xx$x
   y <- xx$y
   r <- xx$coefficients[1]
   a <- xx$coefficients[2]
   b <- xx$coefficients[3]
   if (is.null(xlim)) {    
       xlim <- range(c(a+r*cos(seq(0., 2. * pi, length = 1000.)), x))
   }
   if (is.null(ylim)) {
       ylim <- range(c(b+r*sin(seq(0., 2. * pi, length = 1000.)), y))
   }
   if (is.null(xlab)) xlab <- "X"
   if (is.null(ylab)) ylab <- "Y"
   if (is.null(main)) main <- "Least Square Circle Fit"
   
   midx <- 0.5 * (xlim[2] + xlim[1])
   xlim <- midx + (1 + tol) * 0.5 * c(-1, 1) * (xlim[2] - xlim[1])
   midy <- 0.5 * (ylim[2] + ylim[1])
   ylim <- midy + (1 + tol) * 0.5 * c(-1, 1) * (ylim[2] - ylim[1])
   oldpin <- par("pin")
   xuin <- oxuin <- oldpin[1]/diff(xlim)
   yuin <- oyuin <- oldpin[2]/diff(ylim)
   if (missing(uin)) {
       if (yuin > xuin) xuin <- yuin
       else yuin <- xuin
   } else {
       if (length(uin) == 1) uin <- uin * c(1, 1)
       if (any(c(xuin, yuin) < uin)) stop("uin is too large to fit plot in")
       xuin <- uin[1]; yuin <- uin[2]
   }    
   xlim <- midx + oxuin/xuin * c(-1, 1) * diff(xlim) * 0.5
   ylim <- midy + oyuin/yuin * c(-1, 1) * diff(ylim) * 0.5

   if (!add) {
       plot.default(0, 0, main=main, type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
   }
       lines.default(a+ r*cos(seq(0., 2. * pi, length = 1000.)), b+r*sin(seq(0., 2. * pi, length = 1000.)), ...)
   points.default(x, y, ...)
   text.default(a, b, "+", cex=plus.cex)
   invisible(xx)
}

