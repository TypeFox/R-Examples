# Author: Babak Naimi, naimi.b@gmail.com
# Date :  July 2012
# Version 1.0
# Licence GPL v3

if (!isGeneric("rts")) {
  setGeneric("rts", function(x, time)
    standardGeneric("rts"))
}  

setMethod('rts', signature(x='character'), 
          function(x, time){
            if (!(class(time) %in% c("POSIXct","POSIXt","Date","yearmon","yearqtr"))[1])
              stop("time must be one of time-based classes including: 'POSIXct', 'POSIXt', 'Date', 'yearmon','yearqtr'")
            
            x <- stack(x)
            
            if (nlayers(x) != length(time))
              stop("time must has the same length as the number of layers in RasterStack")
            
            time=xts(1:nlayers(x),time)
            o <- new("RasterStackTS")
            o@raster <- x
            o@time <- time
            o
          }
          )

setMethod('rts', signature(x='RasterStack'), 
          function(x, time){
            if (!(class(time) %in% c("POSIXct","POSIXt","Date","timeDate","yearmon","yearqtr"))[1])
              stop("time must be one of time-based classes including: 'POSIXct', 'POSIXt', 'Date', 'timeDate', 'yearmon','yearqtr'")
            
            if (nlayers(x) != length(time))
              stop("time must has the same length as the number of layers in RasterStack")
            
            time <- xts(1:nlayers(x),time)
            o <- new("RasterStackTS")
            o@raster <- x
            o@time <- time
            o
          }
          )


setMethod('rts', signature(x='RasterBrick'), 
          function(x, time){
            if (!(class(time) %in% c("POSIXct","POSIXt","Date","timeDate","yearmon","yearqtr"))[1])
              stop("time must be one of time-based classes including: 'POSIXct', 'POSIXt', 'Date', 'timeDate', 'yearmon','yearqtr'")
            
            if (nlayers(x) != length(time))
              stop("time must has the same length as the number of layers in RasterBrick")
            
            time=xts(1:nlayers(x),time)
            o <- new("RasterBrickTS")
            o@raster <- x
            o@time <- time
            o
          }
          )
setMethod('rts', signature(x='xts'), 
          function(x){
            structure(x,class=c('rts','xts','zoo'))
          }
          )

setMethod('rts', signature(x='character',time='missing'), 
          function(x){
            if (length(x) > 1) warning("input filename is a vector with more than 1 elements, only the first element is considered!")
            x <- x[1]
            read.rts(x)
          }
          )
