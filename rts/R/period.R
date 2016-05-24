# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2012
# Version 1.0
# Licence GPL v3


if (!isGeneric("period.apply")) {
  setGeneric("period.apply", function(x, INDEX, FUN, ...)
    standardGeneric("period.apply"))
}


setMethod("period.apply", "RasterStackTS",
          function(x,INDEX,FUN,...) {
            FUN <- match.fun(FUN)
            if (length(FUN(1:10)) > 1) stop("Function (FUN) should return one value!") 
            w <- which(!(INDEX %in% 1:length(x@time)))
            if (length(w) > 0) INDEX <- INDEX[-w]
            if (length(INDEX) == 0) stop("INDEX of endpoints is out of range!")
            ep <- sort(INDEX)
            if(ep[length(ep)] != NROW(x@time)) ep <- c(ep,NROW(x@time))
            
            for (i in 1:length(ep)) {
              if (i == 1) {
                if (ep[1] == 1) epr <- 1
                else epr <- 1:ep[1]
                if (length(epr) == 1) xx <- subset(x@raster,1)
                else xx <- calc(subset(x@raster,as.vector(x@time[epr])),FUN, ...)
                xxx <- stack(xx)
                ind <- as.character(index(x@time[ep[1]]))
              } else {
                epr <- (ep[i-1]+1):ep[i]
                if (length(epr) == 1) xx <- subset(x@raster,ep[i])
                else xx <- calc(subset(x@raster,as.vector(x@time[epr])),FUN, ...)
                xxx <- addLayer(xxx,xx)
                ind <- c(ind,as.character(index(x@time[ep[i]])))
              }
            }
            xxx <- rts(xxx,as.POSIXct(ind))
            xxx
          })

setMethod("period.apply", "RasterBrickTS",
          function(x,INDEX,FUN,...) {
            FUN <- match.fun(FUN)
            if (length(FUN(1:10)) > 1) stop("Defined function (FUN) returns more than one value!") 
            w <- which(!(INDEX %in% 1:length(x@time)))
            if (length(w) > 0) INDEX <- INDEX[-w]
            if (length(INDEX) == 0) stop("INDEX of endpoints is out of range!")
            ep <- sort(INDEX)
            if(ep[length(ep)] != NROW(x@time)) ep <- c(ep,NROW(x@time))
            
            for (i in 1:length(ep)) {
              if (i == 1) {
                if (ep[1] == 1) epr <- 1
                else epr <- 1:ep[1]
                if (length(epr) == 1) xx <- subset(x@raster,1)
                else xx <- calc(subset(x@raster,as.vector(x@time[epr])),FUN, ...)
                xxx <- stack(xx)
                ind <- as.character(index(x@time[ep[1]]))
              } else {
                epr <- (ep[i-1]+1):ep[i]
                if (length(epr) == 1) xx <- subset(x@raster,ep[i])
                else xx <- calc(subset(x@raster,as.vector(x@time[epr])),FUN, ...)
                xxx <- addLayer(xxx,xx)
                ind <- c(ind,as.character(index(x@time[ep[i]])))
              }
            }
            xxx <- brick(xxx)
            xxx <- rts(xxx,as.POSIXct(ind))
            xxx
          })