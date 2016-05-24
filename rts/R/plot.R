# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July. 2015
# Version 1.1
# Licence GPL v3


if (!isGeneric("plot")) {
  setGeneric("plot", function(x,y,...)
    standardGeneric("plot"))
}  



setMethod("plot", signature(x='rts',y='ANY'),
          function(x,y,col,main,ylim,...) {
            # part of this function is copied from plot.xts in xts package
            if(nrow(x) < 2) stop("Number of observations should be greater than 1")
            if (missing(y)) y <- 1
            if (y[1] == 'all') y <- 1:ncol(x)
            if (missing(main)) main <- ""
            if (missing(col)) {
              if (length(y) > 1) col <- 1:length(y)
              else col <- 1
            } else if (length(col) != length(y)) col <- rep(col[1],length(y))
            if (missing(ylim)) ylim <- c(min(x[,y],na.rm=TRUE),max(x[,y],na.rm=TRUE))
            
            ep <- axTicksByTime(x, "auto", format.labels = TRUE)
            xycoords <- xy.coords(index(x), x[, y[1]])
            plot(xycoords$x, xycoords$y, type = 'l', axes = FALSE, ann = FALSE, col=col[1], ylim=ylim, ...)
            abline(v = xycoords$x[ep], col = "grey", lty = 4)
            grid(NA, NULL)
            axis(1, at = xycoords$x, labels = FALSE, col = "#BBBBBB")
            axis(1, at = xycoords$x[ep], labels = names(ep), las = 1,lwd = 1, mgp = c(3, 2, 0))
            axis(2)
            box()
            if (length(y) > 1) for (w in 2:length(y)) lines(xycoords$x,x[,y[w]],col=col[w],...)
            title(main)
          }
)


setMethod("plot", signature(x='RasterStackBrickTS','ANY'),
          function(x, y,...) {
            if (missing(y)) y <- as.vector(x@time)
            if (!inherits(try(i <- x@time[y],TRUE), "try-error")) {
              if (length(i) > 0) {
                y <- as.vector(i)
              } else {
                warning("Subscript out of bounds,y is ignored!")
                y <- as.vector(x@time)
              }
            } else {
              warning("No raster is returned for specified time range, y is ignored!")
              y <- as.vector(x@time)
            }
            n <- as.character(index(x@time))[y]
            plot(x=x@raster,y=y,main=n,...)
          }
)
