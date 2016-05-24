#############################################################
#                                                           #
#   points.circular function                                #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: October, 19, 2009                                 #
#   Version: 0.3                                            #
#                                                           #
#   Copyright (C) 2009 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
points.circular <- function(x, pch = 16, cex = 1, stack = FALSE, sep = 0.025, shrink=1, bins=NULL, col=NULL, next.points=NULL, plot.info=NULL, zero=NULL, rotation=NULL, ...) {
   if (is.matrix(x) | is.data.frame(x)) {
      nseries <- ncol(x)
   } else {
      nseries <- 1
   }
   xx <- as.data.frame(x)
  
   xcircularp <- attr(as.circular(xx[,1]), "circularp")
   type <- xcircularp$type
   modulo <- xcircularp$modulo
   if (is.null(plot.info)) {
      if (is.null(zero))
         zero <- xcircularp$zero
      if (is.null(rotation))
         rotation <- xcircularp$rotation
      if (is.null(next.points))
         next.points <- 0
   } else {
      zero <- plot.info$zero
      rotation <- plot.info$rotation
      if (is.null(next.points))
         next.points <- plot.info$next.points
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
            
   for (iseries in 1:nseries) {
      x <- xx[,iseries]
      x <- na.omit(x)
      n <- length(x)
      if (n) {
         x <- conversion.circular(x, units="radians")
         attr(x, "circularp") <- attr(x, "class") <- NULL
         if (rotation=="clock")
            x <- -x
         x <- x+zero
         x <- x%%(2*pi)
         PointsCircularRad(x, bins, stack, col, pch, iseries, nseries, sep, next.points, shrink, cex, ...) 
      }
   }
return(invisible(list(zero=zero, rotation=rotation, next.points=next.points+nseries*sep)))
}

PointsCircularRad <- function(x, bins, stack, col, pch, iseries, nseries, sep, next.points, shrink, cex, ...) {
#### x musts be in modulo 2pi  
   if (!stack) {
      z <- cos(x)
      y <- sin(x)
      r <- 1+((iseries-1)*sep+next.points)*shrink
      points.default(z*r, y*r, cex=cex, pch=pch[iseries], col = col[iseries], ...)
   } else {
      x[x >= 2*pi] <- 2*pi-4*.Machine$double.eps
      arc <- (2 * pi)/bins
      pos.bins <- ((1:nseries)-1/2)*arc/nseries-arc/2
#      bins.count <- c(1:bins)
#      for (i in 1:bins) {
#         bins.count[i] <- sum(x < i * arc & x >= (i - 1) * arc)
#      }
      breaks <- seq(0,2*pi,length.out=(bins+1))
      bins.count <- hist.default(x, breaks=breaks, plot=FALSE, right=TRUE)$counts
###### TO BE USED IN THE FUTURE .C("bincount", x, as.integer(length(x)), seq(0,2*pi,length.out=bins), as.integer(bins+1), counts = integer(bins), right = as.logical(TRUE), include = as.logical(FALSE), naok = FALSE, NAOK = FALSE, DUP = FALSE, PACKAGE = "base")$counts
      mids <- seq(arc/2, 2 * pi - pi/bins, length = bins) + pos.bins[iseries]
      index <- cex*sep
      for (i in 1:bins) {
         if (bins.count[i] != 0) {
            for (j in 0:(bins.count[i] - 1)) {
               r <- 1 + j * index
               z <- r * cos(mids[i])
               y <- r * sin(mids[i])
               points.default(z, y, cex=cex, pch=pch[iseries], col=col[iseries], ...)
            }
         }
      }
   }
}

