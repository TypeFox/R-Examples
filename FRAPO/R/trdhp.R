##
##
## Generic for trdhp trend 
##
setGeneric(name = "trdhp", def = function(y, lambda){standardGeneric("trdhp")})
##
## Methods for trdhp
## (in alphabetical order)
##
## for class data.frame
setMethod(f = "trdhp",
          signature = c(y = "data.frame"),
          definition = function(y, lambda){
            if(!all(apply(y, 2, is.numeric))){
              stop("Data frame has non-numeric columns.\n")
            }
            trd <- apply(y, 2, trdhp, lambda = lambda)
            return(trd)
          }
)
## for class matrix
setMethod(f = "trdhp",
          signature = c(y = "matrix"),
          definition = function(y, lambda){
            trd <- apply(y, 2, trdhp, lambda = lambda)
            return(trd)
          }
)
## for class mts
setMethod(f = "trdhp",
          signature = c(y = "mts"),
          definition = function(y, lambda){
            trd <- apply(y, 2, trdhp, lambda = lambda)
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
## for class numeric
setMethod(f = "trdhp",
          signature = c(y = "numeric"),
          definition = function(y, lambda){
            ident <- diag(length(y))
            d <- diff(ident, d = 2)
            trd <- solve(ident + lambda * crossprod(d), y)
            return(trd)
          }
)
## for class timeSeries
setMethod(f = "trdhp",
          signature = c(y = "timeSeries"),
          definition = function(y, lambda){
            trd <- apply(y, 2, trdhp, lambda = lambda)
            return(trd)
          }
)
## for class ts
setMethod(f = "trdhp",
          signature = c(y = "ts"),
          definition = function(y, lambda){
            trd <- trdhp(c(y), lambda = lambda)
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
## for class xts
setMethod(f = "trdhp",
          signature = c(y = "xts"),
          definition = function(y, lambda){
            trd <- trdhp(coredata(y), lambda = lambda)
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
## for class zoo
setMethod(f = "trdhp",
          signature = c(y = "zoo"),
          definition = function(y, lambda){
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
