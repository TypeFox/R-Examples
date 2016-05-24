# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2015
# Version 1.1
# Licence GPL v3

setMethod("[[", c("RasterStackBrickTS","ANY","ANY"),
          function(x,i,drop=TRUE, ...) {
            if ( missing(i)) { stop('you must provide an index') }
            if (!inherits(try(i <- x@time[i],T), "try-error")) {
              if (length(i) > 1) {
                rts(subset(x@raster,as.vector(i),drop=drop,...),index(i))
              } else {
                if (length(i) == 1) {
                  x <- subset(x@raster,as.vector(i),drop=drop,...)
                  names(x) <- as.character(index(i))
                  x
                } else stop("There is no data for specified time range!")
              }
            } else {
              stop("There is no data for specified time range, or subscript out of bounds")
            }
          })

if (!isGeneric("subset")) {
  setGeneric("subset", function(x, ...)
    standardGeneric("subset"))
}

setMethod("subset","RasterStackBrickTS",
          function(x, subset, drop=TRUE, ...) {
            if ( missing(subset)) { stop('you must provide an index') }
            x[[subset,drop=drop,...]]
          })
