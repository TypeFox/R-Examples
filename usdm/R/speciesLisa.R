# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Sep 2012
# Version 1.0
# Licence GPL v3


if (!isGeneric("speciesLisa")) {
  setGeneric("speciesLisa", function(x, y, uncertainty, statistic="K1",weights)
    standardGeneric("speciesLisa"))
}  
setMethod('speciesLisa', signature(x='Raster',y='SpatialPoints'), 
          function(x, y, uncertainty, statistic="K1",weights) {
            if (nlayers(x) > 1) {
              if (missing(weights)) weights <- rep(1/nlayers(x),nlayers(x))
              else if (length(weights) != nlayers(x)) stop("the length of weights should be equal to the number of layers in the raster object")
              
              weights <- weights/sum(weights)
              o <- lisa(x=x,y=y,d1=0,d2=uncertainty,statistic=statistic)
              n <- new("speciesLISA")
              n@species <- y
              n@data <- NA
              n@LISAs <- o
              n@weights <- weights
              n@statistic <- statistic
              for (i in 1:length(weights)) o[,i] <- o[,i] * weights[i]
              n@LISA <- apply(o,1,sum)
            } else {
              o <- lisa(x=x,y=y,d1=0,d2=uncertainty,statistic=statistic)
              n <- new("speciesLISA")
              n@species <- y
              n@LISAs <- o
              n@weights <- NA
              n@statistic <- statistic
              n@LISA <- o
            }
            n
          }
)

setMethod('speciesLisa', signature(x='Raster',y='SpatialPointsDataFrame'), 
          function(x, y, uncertainty, statistic="K1",weights) {
            if (nlayers(x) > 1) {
              if (missing(weights)) weights <- rep(1/nlayers(x),nlayers(x))
              else if (length(weights) != nlayers(x)) stop("the length of weights should be equal to the number of layers in the raster object")
              
              weights <- weights/sum(weights)
              o <- lisa(x=x,y=y,d1=0,d2=uncertainty,statistic=statistic)
              n <- new("speciesLISA")
              n@species <- as(y,"SpatialPoints")
              n@data <- as(y,'data.frame')
              n@LISAs <- o
              n@weights <- weights
              n@statistic <- statistic
              for (i in 1:length(weights)) o[,i] <- o[,i] * weights[i]
              n@LISA <- apply(o,1,sum)
            } else {
              o <- lisa(x=x,y=y,d1=0,d2=uncertainty,statistic=statistic)
              n <- new("speciesLISA")
              n@species <- as(y,"SpatialPoints")
              n@data <- as(y,'data.frame')
              n@LISAs <- o
              n@statistic <- statistic
              n@weights <- NA
              n@LISA <- o
            }
            n
          }
)