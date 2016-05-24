##
##
## Generic for trdwma trend 
##
setGeneric(name = "trdwma", def = function(y, weights, trim = FALSE){standardGeneric("trdwma")})
##
## Methods for trdwma
## (in alphabetical order)
##
## for class data.frame
setMethod(f = "trdwma",
          signature = c(y = "data.frame"),
          definition = function(y,  weights, trim = FALSE){
            if(!all(apply(y, 2, is.numeric))){
              stop("Data frame has non-numeric columns.\n")
            }
            trd <- as.data.frame(apply(y, 2, trdwma, weights = weights, trim = trim), ncol = ncol(y))
            colnames(trd) <- colnames(y)
            ifelse(trim, rownames(trd) <- rownames(y)[-c(1:(weights - 1))], rownames(trd) <- rownames(y))             
            return(as.data.frame(trd))
          }
)
## for class matrix
setMethod(f = "trdwma",
          signature = c(y = "matrix"),
          definition = function(y, weights, trim = FALSE){
            trd <- matrix(apply(y, 2, trdwma, weights = weights, trim = trim), ncol = ncol(y))
            colnames(trd) <- colnames(y)
            return(trd)
          }
)
## for class mts
setMethod(f = "trdwma",
          signature = c(y = "mts"),
          definition = function(y, weights, trim = FALSE){
            trd <- matrix(apply(y, 2, trdwma, weights = weights, trim = FALSE), ncol = ncol(y))
            attributes(trd) <- attributes(y)
            if(trim) trd <- window(trd, start = time(y)[length(weights)])
            return(trd)
          }
)
## for class numeric
setMethod(f = "trdwma",
          signature = c(y = "numeric"),
          definition = function(y, weights, trim = FALSE){
            weights <- as.numeric(weights)
            if(!(sum(weights) == 1)){
              warning("\nThe sum of the weights is not equal to one.\n")            
            }
            if(length(weights) > length(y)){
              stop("\nNumber of weights is greater than length of series.\n")
            }
            trd <- c(filter(y, filter = weights, sides = 1))              
            if(trim){
              trd <- trd[-c(1:(length(weights) -1))]
            }
            return(trd)
          }
)
## for class timeSeries
setMethod(f = "trdwma",
          signature = c(y = "timeSeries"),
          definition = function(y, weights, trim = FALSE){
            trd <- apply(y, 2, trdwma, weights = weights, trim = FALSE)
            if(trim) trd <- window(trd, start = time(y)[length(weights)], end = end(y)) 
            return(trd)
          }
)
## for class ts
setMethod(f = "trdwma",
          signature = c(y = "ts"),
          definition = function(y, weights, trim = FALSE){
            trd <- trdwma(c(y), weights = weights, trim = FALSE)
            attributes(trd) <- attributes(y)
            if(trim) trd <- window(trd, start = time(y)[length(weights)])
            return(trd)
          }
)
## for class xts
setMethod(f = "trdwma",
          signature = c(y = "xts"),
          definition = function(y, weights, trim = FALSE){
            yc <- as.matrix(coredata(y))
            trd <- trdwma(yc, weights = weights, trim = FALSE)
            attributes(trd) <- attributes(y)
            if(trim) trd <- window(trd, start = index(y)[length(weights)])
            return(trd)
          }
)
## for class zoo
setMethod(f = "trdwma",
          signature = c(y = "zoo"),
          definition = function(y, weights, trim = FALSE){
            yc <- as.matrix(coredata(y))
            trd <- trdwma(yc, weights = weights, trim = FALSE)
            attributes(trd) <- attributes(y)
            if(trim) trd <- window(trd, start = index(y)[length(weights)])
            return(trd)
          }
)
