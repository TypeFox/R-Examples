
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
#   Date: March, 15, 2011                                   #
#   Copyright (C) 2011 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-3                                           #
#                                                           #
#   Modified by Hiroyoshi Arai                              #
#   Date: October, 19, 2010                                 #
#   Added arguments                                         #
#     upper: if TRUE, upper-closed (lower-open) intervals.  #
#     radii.scale: "sqrt"(default) or "linear"              #
#     border: the color to draw the border.                 #
#   Revised argument                                        #
#     col: the color to filling the sector.                 #
#                                                           #
#############################################################

rose.diag <- function(x, pch = 16, cex=1, axes = TRUE, shrink = 1, bins=NULL, upper=TRUE, ticks = TRUE, tcl=0.025, tcl.text=0.125, radii.scale = c("sqrt", "linear"), border=NULL, col=NULL, tol = 0.04, uin=NULL, xlim=c(-1, 1), ylim=c(-1, 1), prop = 1, digits=2, plot.info=NULL, units=NULL, template=NULL, zero=NULL, rotation=NULL, main=NULL, sub=NULL, xlab="", ylab="", add=FALSE, control.circle = circle.control(), ...) {

   radii.scale <- match.arg(radii.scale)
   
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
      CirclePlotRad(xlim=xlim, ylim=ylim, uin=uin, shrink=shrink, tol=tol, main=main, sub=sub, xlab=xlab, ylab=ylab, control.circle=control.circle)
   }
   
   if (is.null(bins)) {
      bins <- NROW(x)
   } else {
      bins <- round(bins)
      if (bins<=0)
         stop("bins must be non negative")
   }

   if (is.null(border)) {
      border <- seq(nseries)
   } else {
      if (length(border)!=nseries) {
         border <- rep(border, nseries)[1:nseries]
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
        # x <- x+zero
         if (template=="clock12")
           x <- 2*x
         x <- x%%(2*pi)
         RosediagRad(x, zero=zero, rotation, bins, upper, radii.scale, prop, border[iseries], col, ...)
      }
   }
   return(invisible(list(zero=zero, rotation=rotation, next.points=0)))    
}

RosediagRad <- function(x, zero, rotation, bins, upper, radii.scale, prop, border, col, ...) {
#### x musts be in modulo 2pi
    n <- length(x)
    freq <- rep(0, bins)
    arc <- (2 * pi)/bins
    if (!is.logical(upper))
       stop("upper must be logical")
    if (upper == TRUE)
       x[x == 0] <- 2*pi
    
    x[x >= 2*pi] <- 2*pi-4*.Machine$double.eps
#    for (i in 1:bins) {
#       freq[i] <- sum(x < i * arc & x >= (i - 1) * arc)
#    }
    breaks <- seq(0,2*pi,length.out=(bins+1))
    freq <- hist.default(x, breaks=breaks, plot=FALSE, right=upper)$counts   
    rel.freq <- freq/n
    if (rotation == "clock")
       rel.freq <- rev(rel.freq)
    
    if (radii.scale == "sqrt") {
       radius <- sqrt(rel.freq)*prop
    } else {
       radius <- rel.freq*prop
    }
    sector <- seq(0, 2 * pi - (2 * pi)/bins, length = bins)
    mids <- seq(arc/2, 2 * pi - pi/bins, length = bins)
    for (i in 1:bins) {
       if (rel.freq[i] != 0) {
          xx <- c(0, radius[i]*cos(seq(sector[i], sector[i]+(2*pi)/bins, length=1000/bins)+zero), 0)
          yy <- c(0, radius[i]*sin(seq(sector[i], sector[i]+(2*pi)/bins, length=1000/bins)+zero), 0)
          polygon(xx, yy, border=border, col=col, ...)
       }
    }
}
