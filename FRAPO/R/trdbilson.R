##
##
## Generic for trdbilson trend 
##
setGeneric(name = "trdbilson", def = function(y, exponent){standardGeneric("trdbilson")})
##
## Methods for trdbilson
## (in alphabetical order)
##
## for class data.frame
setMethod(f = "trdbilson",
          signature = c(y = "data.frame"),
          definition = function(y, exponent){
            if(!all(apply(y, 2, is.numeric))){
              stop("Data frame has non-numeric columns.\n")
            }
            trd <- apply(y, 2, trdbilson, exponent = exponent)
            return(trd)
          }
)
## for class matrix
setMethod(f = "trdbilson",
          signature = c(y = "matrix"),
          definition = function(y, exponent){
            trd <- apply(y, 2, trdbilson, exponent = exponent)
            return(trd)
          }
)
## for class mts
setMethod(f = "trdbilson",
          signature = c(y = "mts"),
          definition = function(y, exponent){
            trd <- apply(y, 2, trdbilson, exponent)
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
## for class numeric
setMethod(f = "trdbilson",
          signature = c(y = "numeric"),
          definition = function(y, exponent){
            absval <- abs(y)
            sigval <- sign(y)
            trd <- sigval * absval^(1 - absval^exponent)
            return(trd)
          }
)
## for class timeSeries
setMethod(f = "trdbilson",
          signature = c(y = "timeSeries"),
          definition = function(y, exponent){
            trd <- apply(y, 2, trdbilson, exponent)
            return(trd)
          }
)
## for class ts
setMethod(f = "trdbilson",
          signature = c(y = "ts"),
          definition = function(y, exponent){
            trd <- trdbilson(c(y), exponent)
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
## for class xts
setMethod(f = "trdbilson",
          signature = c(y = "xts"),
          definition = function(y, exponent){
            yc <- as.matrix(coredata(y))
            trd <- apply(yc, 2, trdbilson, exponent = exponent)
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
## for class zoo
setMethod(f = "trdbilson",
          signature = c(y = "zoo"),
          definition = function(y, exponent){
            yc <- as.matrix(coredata(y))
            trd <- apply(yc, 2, trdbilson, exponent = exponent)            
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
