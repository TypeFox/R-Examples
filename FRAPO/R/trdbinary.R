##
##
## Generic for trdbinary trend 
##
setGeneric(name = "trdbinary", def = function(y){standardGeneric("trdbinary")})
##
## Methods for trdbinary
## (in alphabetical order)
##
## for class data.frame
setMethod(f = "trdbinary",
          signature = c(y = "data.frame"),
          definition = function(y){
            if(!all(apply(y, 2, is.numeric))){
              stop("Data frame has non-numeric columns.\n")
            }
            trd <- apply(y, 2, trdbinary)
            return(trd)
          }
)
## for class matrix
setMethod(f = "trdbinary",
          signature = c(y = "matrix"),
          definition = function(y){
            trd <- apply(y, 2, trdbinary)
            return(trd)
          }
)
## for class mts
setMethod(f = "trdbinary",
          signature = c(y = "mts"),
          definition = function(y){
            trd <- apply(y, 2, trdbinary)
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
## for class numeric
setMethod(f = "trdbinary",
          signature = c(y = "numeric"),
          definition = function(y){
            trd <- sign(y) * pmin(abs(4 / pi * atan(y)), 1)
            return(trd)
          }
)
## for class timeSeries
setMethod(f = "trdbinary",
          signature = c(y = "timeSeries"),
          definition = function(y){
            trd <- apply(y, 2, trdbinary)
            return(trd)
          }
)
## for class ts
setMethod(f = "trdbinary",
          signature = c(y = "ts"),
          definition = function(y){
            trd <- trdbinary(c(y))
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
## for class xts
setMethod(f = "trdbinary",
          signature = c(y = "xts"),
          definition = function(y){
            yc <- as.matrix(coredata(y))
            trd <- apply(yc, 2, trdbinary)
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
## for class zoo
setMethod(f = "trdbinary",
          signature = c(y = "zoo"),
          definition = function(y){
            yc <- as.matrix(coredata(y))
            trd <- apply(yc, 2, trdbinary)            
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
