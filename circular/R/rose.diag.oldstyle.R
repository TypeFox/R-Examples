#####This code is not used anymore. It is here for historical reason. Please use
#####the code in the file rose.diag.R

#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   rose.diag function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: October, 18, 2009                                 #
#   Copyright (C) 2009 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-2                                           #
#                                                           #
#############################################################

rose.diag.oldstyle <- function(x, pch = 16, cex=1, axes = TRUE, shrink = 1, bins=NULL, ticks = TRUE, tcl=0.025, tcl.text=0.125, col=NULL, tol = 0.04, uin=NULL, xlim=c(-1, 1), ylim=c(-1, 1), prop = 1, digits=2, plot.info=NULL, units=NULL, template=NULL, zero=NULL, rotation=NULL, main="", xlab="", ylab="", add=FALSE, ...) {
  
   if (is.matrix(x) | is.data.frame(x)) {
      nseries <- ncol(x)
   } else {
      nseries <- 1
   }
   xx <- as.data.frame(x)
  
   xcircularp <- attr(as.circular(xx[,1]), "circularp")
#   type <- xcircularp$type
   modulo <- xcircularp$modulo
   if (is.null(units)) 
      units <- xcircularp$units
   if (is.null(plot.info)) {
      if (is.null(template))
         template <- xcircularp$template
      if (template=="geographics" | template=="clock24") {
         zero <- pi/2
         rotation <- "clock"
      } else if (template=="clock12") {
         zero <- pi/2
         rotation <- "clock"
         modulo <- "pi"
      } else {
         if (is.null(zero))
            zero <- xcircularp$zero
         if (is.null(rotation))
            rotation <- xcircularp$rotation
      }
      next.points <- 0
   } else {
      zero <- plot.info$zero
      rotation <- plot.info$rotation
      next.points <- plot.info$next.points
   }
   
   if (!add) {
      CirclePlotRad(xlim, ylim, uin, shrink, tol, 1000, main=main, xlab=xlab, ylab=ylab)
   }
   
   if (is.null(bins)) {
      bins <- NROW(x)
   } else {
      bins <- round(bins)
      if (bins<=0)
         stop("bins must be non negative")
   }

   if (is.null(col)) {
      col <- seq(nseries)
   } else {
      if (length(col)!=nseries) {
         col <- rep(col, nseries)[1:nseries]
      }
   }
   pch <- rep(pch, nseries, length.out=nseries)

   if (axes) {
      axis.circular(units=units, template=template, zero=zero, rotation=rotation, digits=digits, cex=cex, tcl=tcl, tcl.text=tcl.text)
   }

   if (!is.logical(ticks))
      stop("ticks must be logical")
      
   if (ticks) {
      at <- circular((0:bins)/bins*2*pi, zero=zero, rotation=rotation)
      ticks.circular(at, tcl=tcl)
   }
      
   for (iseries in 1:nseries) {   
      x <- xx[,iseries]
      x <- na.omit(x)
      n <- length(x)      
      if (n) {
         x <- conversion.circular(x, units="radians", modulo=modulo)
         attr(x, "circularp") <- attr(x, "class") <- NULL
         if (rotation=="clock")
            x <- -x
         x <- x+zero
         if (template=="clock12")
           x <- 2*x
         x <- x%%(2*pi)
         RosediagOSRad(x, bins, prop, col[iseries], ...)
      }
   }
   return(invisible(list(zero=zero, rotation=rotation, next.points=0)))    
}

RosediagOSRad <- function(x, bins, prop, col, ...) {
#### x musts be in modulo 2pi
    n <- length(x)
    freq <- rep(0, bins)
    arc <- (2 * pi)/bins
    x[x >= 2*pi] <- 2*pi-4*.Machine$double.eps
#    for (i in 1:bins) {
#       freq[i] <- sum(x < i * arc & x >= (i - 1) * arc)
#    }
    breaks <- seq(0,2*pi,length.out=(bins+1))
    freq <- hist.default(x, breaks=breaks, plot=FALSE, right=TRUE)$counts   
    rel.freq <- freq/n
    radius <- sqrt(rel.freq) * prop
    sector <- seq(0, 2 * pi - (2 * pi)/bins, length = bins)
    mids <- seq(arc/2, 2 * pi - pi/bins, length = bins)
    for (i in 1:bins) {
       if (rel.freq[i] != 0) {
          lines.default(c(0, radius[i] * cos(sector[i])), c(0, radius[i] * sin(sector[i])), col=col, ...)
          lines.default(c(0, radius[i] * cos(sector[i] + (2 * pi)/bins)), c(0, radius[i] * sin(sector[i] + (2 * pi)/bins)), col=col, ...)
          lines.default(c(radius[i] * cos(sector[i]), radius[i] * cos(sector[i] + (2 * pi)/bins)), c(radius[i] * sin(sector[i]), radius[i] * sin(sector[i] + (2 * pi)/bins)), col=col, ...)
       }
    }
}
