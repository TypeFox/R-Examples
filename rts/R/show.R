# Author: Babak Naimi, naimi.b@gmail.com
# This is based on 'show' function from raster package 
# Date :  November 2012
# Version 1.1
# Licence GPL v3

setMethod ('show' , 'RasterStackBrickTS',
           function ( object ) {
             if (length(object@time) > 1) {
               p <- periodicity(object@time)
               cat ('Raster Time Series with',p$scale, 'periodicity from',as.character(p$start),'to',as.character(p$end),'\n')
             } 
             cat ('class       :' , class ( object ) , '\n')
                          
             if (filename(object@raster) != '') {
               cat ('raster filename    :' , filename(object@raster), '\n')
             }
             nl <- nlayers(object@raster)
             if (nl == 0) {
               cat ('nlayers     :' , nl, '\n')
             } else {
               cat ('raster dimensions  : ', nrow(object@raster), ', ', ncol(object@raster), ', ', ncell(object@raster), ', ', nl, '  (nrow, ncol, ncell, nlayers)\n', sep="" ) 
               cat ('raster resolution  : ' , xres(object@raster), ', ', yres(object@raster), '  (x, y)\n', sep="")
               cat ('raster extent      : ' , object@raster@extent@xmin, ', ', object@raster@extent@xmax, ', ', object@raster@extent@ymin, ', ', object@raster@extent@ymax, '  (xmin, xmax, ymin, ymax)\n', sep="")
               cat ('coord. ref. :' , projection(object@raster, TRUE), '\n')
               
               minv <- format(minValue(object@raster), digits=2)
               maxv <- format(maxValue(object@raster), digits=2)
               minv <- gsub('Inf', '?', minv)
               maxv <- gsub('-Inf', '?', maxv)
               if (nl > 10) {
                 minv <- c(minv[1:10], '...')
                 maxv <- c(maxv[1:10], '...')
               }
               cat('min values  :', paste(minv, collapse=' '), '\n')
               cat('max values  :', paste(maxv, collapse=' '), '\n')
             }
             
             z <- getZ(object@raster)
             if (length(z) > 0) {
               name <- object@raster@zname
               if (name == '') name <- 'z-value'
               name <- paste(sprintf("%-12s", name), ':', sep='')
               if (length(z) < 10) {
                 cat(name, paste(z, collapse=', '), '\n')
               } else {
                 z <- summary(z)
                 cat(name, paste(z, collapse=' ... '), '(summary)\n')
               }
             }
             
             cat ('\n')
           }
           )
