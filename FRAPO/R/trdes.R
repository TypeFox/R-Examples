##
##
## Generic for trdes trend 
##
setGeneric(name = "trdes", def = function(y, lambda, init = NULL){standardGeneric("trdes")})
##
## Methods for trdes
## (in alphabetical order)
##
## for class data.frame
setMethod(f = "trdes",
          signature = c(y = "data.frame"),
          definition = function(y, lambda, init = NULL){
            if(!all(apply(y, 2, is.numeric))){
              stop("Data frame has non-numeric columns.\n")
            }
            if(is.null(init)){
              trd <- apply(y, 2, trdes, lambda = lambda)
            } else {
              idx <- 1:ncol(y)
              init <- as.vector(init)
              trd <- sapply(idx, function(i) trdes(y[, i], lambda = lambda, init = init[i]))
            }
            return(trd)
          }
)
## for class matrix
setMethod(f = "trdes",
          signature = c(y = "matrix"),
          definition = function(y, lambda, init = NULL){
            if(is.null(init)){
              trd <- apply(y, 2, trdes, lambda = lambda)
            } else {
              idx <- 1:ncol(y)
              init <- as.vector(init)
              trd <- sapply(idx, function(i) trdes(y[, i], lambda = lambda, init = init[i]))
            }
            return(trd)
          }
)
## for class mts
setMethod(f = "trdes",
          signature = c(y = "mts"),
          definition = function(y, lambda, init = NULL){
            if(is.null(init)){
              trd <- apply(y, 2, trdes, lambda = lambda)
            } else {
              idx <- 1:ncol(y)
              init <- as.vector(init)
              trd <- sapply(idx, function(i) trdes(y[, i], lambda = lambda, init = init[i]))
            }
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
## for class numeric
setMethod(f = "trdes",
          signature = c(y = "numeric"),
          definition = function(y, lambda, init = NULL){
            if((0 > lambda) | (lambda > 1)){
              stop("\nThe parameter lambda must be in the interval (0, 1).\n")
            }
            if(is.null(init)){
              trd <- filter(x = lambda * y, filter = 1 - lambda, method = "recursive")
            } else {
              trd <- filter(x = lambda * y, filter = 1 - lambda, method = "recursive", init = init)
            }
            return(trd)
          }
)
## for class timeSeries
setMethod(f = "trdes",
          signature = c(y = "timeSeries"),
          definition = function(y, lambda, init = NULL){
            if(is.null(init)){
              trd <- apply(y, 2, trdes, lambda = lambda)
            } else {
              idx <- 1:ncol(y)
              init <- as.vector(init)
              trd <- sapply(idx, function(i) trdes(c(y[, i]), lambda = lambda, init = init[i]))
              trd <- timeSeries(trd, charvec = time(y))
              attributes(trd) <- attributes(y)
            }
            return(trd)
          }
)
## for class ts
setMethod(f = "trdes",
          signature = c(y = "ts"),
          definition = function(y, lambda, init = NULL){
            if(is.null(init)){
              trd <- trdes(c(y), lambda = lambda)
            } else {
              trd <- trdes(c(y), lambda = lambda, init = init)
            }
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
## for class xts
setMethod(f = "trdes",
          signature = c(y = "xts"),
          definition = function(y, lambda, init = NULL){
            yc <- as.matrix(coredata(y))
            if(is.null(init)){
              trd <- trdes(yc, lambda = lambda)
            } else {
              idx <- 1:ncol(yc)
              init <- as.vector(init)
              trd <- sapply(idx, function(i) trdes(yc[, i], lambda = lambda, init = init[i]))
            }
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
## for class zoo
setMethod(f = "trdes",
          signature = c(y = "zoo"),
          definition = function(y, lambda, init = NULL){
            yc <- as.matrix(coredata(y))
            if(is.null(init)){
              trd <- trdes(yc, lambda = lambda)
            } else {
              idx <- 1:ncol(yc)
              init <- as.vector(init)
              trd <- sapply(idx, function(i) trdes(yc[, i], lambda = lambda, init = init[i]))
            }
            attributes(trd) <- attributes(y)
            return(trd)
          }
)
