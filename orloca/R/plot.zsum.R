#
# Plots for objective function
#

contour.loca.p <- function(x, lp=numeric(0), xmin=min(min(x@x), xleft), xmax=max(max(x@x), xright), ymin=min(min(x@y), ybottom), ymax=max(max(x@y), ytop), n=100, img=NULL, xleft=min(x@x), ybottom=min(x@y), xright=max(x@x), ytop=max(x@y), ...)
   {
   .x<-seq(xmin, xmax, length.out=n)
   .y<-seq(ymin, ymax, length.out=n)
   .z<-matrix(1, nrow=n, ncol=n, byrow=TRUE)
   if (length(lp) == 0)
     {
       for(i in 1:n)
         for(j in 1:n)
           .z[i,j] <- zsum(x, .x[i], .y[j])
     }
   else if (lp >= 1)
     {
       for(i in 1:n)
         for(j in 1:n)
           .z[i,j] <- zsumlp(x, .x[i], .y[j], p=lp)
    }
   else stop(paste(lp, gettext("is not a valid value for lp, use 1 <= lp", domain = "R-orloca")))
   contour(.x, .y, .z, ...) 
   if (!is.null(img)) {
     if (is.raster(.img <- img) || is.raster(.img <- as.raster(img))) {
       require('png')
       rasterImage(.img, xleft, ybottom, xright, ytop)
       contour(.x, .y, .z, add=TRUE, ...) 
       }
     else warning(gettext("The given img object is not a raster image and cannot be coerce to it.", domain = "R-orloca"))
   }
   invisible(1)
   }

persp.loca.p <- function(x, lp=numeric(0), xmin=min(x@x), xmax=max(x@x), ymin=min(x@y), ymax=max(x@y), n=100, ...)
   {
   .x<-seq(xmin, xmax, length.out=n)
   .y<-seq(ymin, ymax, length.out=n)
   .z<-matrix(1, nrow=n, ncol=n, byrow=TRUE)
   if (length(lp) == 0)
     {
       for(i in 1:n)
         for(j in 1:n)
           .z[i,j] <- zsum(x, .x[i], .y[j])
     }
   else if (lp >= 1)
     {
       for(i in 1:n)
         for(j in 1:n)
           .z[i,j] <- zsumlp(x, .x[i], .y[j], p=lp)
    }
   else stop(paste(lp, gettext("is not a valid value for lp, use 1 <= lp", domain = "R-orloca")))
   persp(.x, .y, .z, ...)
   invisible(1)
   }

