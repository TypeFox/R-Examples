#
# Graphics for loca.p
#
plot.loca.p <- function(x, xlab="", ylab="", main=gettext("Plot of loca.p object", domain = "R-orloca"), img=NULL, xlim=c(min(xleft, min(x@x)), max(xright, max(x@x))), ylim=c(min(ybottom, min(x@y)), max(ytop, max(x@y))), xleft=min(x@x), ybottom=min(x@y), xright=max(x@x), ytop=max(x@y), ...)
   {
   plot(x@x, x@y, xlab=xlab, ylab=ylab, main=main, xlim=xlim, ylim=ylim, ...)
   if (!is.null(img)) {
     if (is.raster(.img <- img) || is.raster(.img <- as.raster(img))) {
       require('png')
       rasterImage(.img, xleft, ybottom, xright, ytop)
       }
     else warning(gettext("The given img object is not a raster image and cannot be coerce to it.", domain = "R-orloca"))
     }
   points(x@x, x@y)
   }
