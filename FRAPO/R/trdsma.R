##
##
## Generic for trdsma trend 
##
setGeneric(name = "trdsma", def = function(y, n.periods, trim = FALSE){standardGeneric("trdsma")})
##
## Methods for trdsma
## (in alphabetical order)
##
## for class data.frame
setMethod(f = "trdsma",
          signature = c(y = "data.frame"),
          definition = function(y,  n.periods, trim = FALSE){
            if(!all(apply(y, 2, is.numeric))){
              stop("Data frame has non-numeric columns.\n")
            }
            trd <- as.data.frame(apply(y, 2, trdsma, n.periods = n.periods, trim = trim), ncol = ncol(y))
            colnames(trd) <- colnames(y)
            ifelse(trim, rownames(trd) <- rownames(y)[-c(1:(n.periods - 1))], rownames(trd) <- rownames(y))             
            return(as.data.frame(trd))
          }
)
## for class matrix
setMethod(f = "trdsma",
          signature = c(y = "matrix"),
          definition = function(y, n.periods, trim = FALSE){
            trd <- matrix(apply(y, 2, trdsma, n.periods = n.periods, trim = trim), ncol = ncol(y))
            colnames(trd) <- colnames(y)
            return(trd)
          }
)
## for class mts
setMethod(f = "trdsma",
          signature = c(y = "mts"),
          definition = function(y, n.periods, trim = FALSE){
            trd <- matrix(apply(y, 2, trdsma, n.periods = n.periods, trim = FALSE), ncol = ncol(y))
            attributes(trd) <- attributes(y)
            if(trim) trd <- window(trd, start = time(y)[n.periods])
            return(trd)
          }
)
## for class numeric
setMethod(f = "trdsma",
          signature = c(y = "numeric"),
          definition = function(y, n.periods, trim = FALSE){
            n.periods <- abs(as.integer(n.periods))
            if(n.periods > length(y)){
              stop("\nNumber of periods is greater than length of series.\n")
            }
            trd <- c(filter(y, filter = rep(1 / n.periods, n.periods), sides = 1))              
            if(trim){
              trd <- na.omit(trd)
            }
            return(trd)
          }
)
## for class timeSeries
setMethod(f = "trdsma",
          signature = c(y = "timeSeries"),
          definition = function(y, n.periods, trim = FALSE){
            trd <- apply(y, 2, trdsma, n.periods = n.periods, trim = FALSE)
            trd <- timeSeries(trd, charvec = time(y))
            if(trim) trd <- window(trd, start = time(y)[n.periods], end = end(y)) 
            return(trd)
          }
)
## for class ts
setMethod(f = "trdsma",
          signature = c(y = "ts"),
          definition = function(y, n.periods, trim = FALSE){
            trd <- trdsma(c(y), n.periods = n.periods, trim = FALSE)
            attributes(trd) <- attributes(y)
            if(trim) trd <- window(trd, start = time(y)[n.periods])
            return(trd)
          }
)
## for class xts
setMethod(f = "trdsma",
          signature = c(y = "xts"),
          definition = function(y, n.periods, trim = FALSE){
            yc <- as.matrix(coredata(y))
            trd <- trdsma(yc, n.periods = n.periods, trim = FALSE)
            attributes(trd) <- attributes(y)
            if(trim) trd <- window(trd, start = index(y)[n.periods])
            return(trd)
          }
)
## for class zoo
setMethod(f = "trdsma",
          signature = c(y = "zoo"),
          definition = function(y, n.periods, trim = FALSE){
            yc <- as.matrix(coredata(y))
            trd <- trdsma(yc, n.periods = n.periods, trim = FALSE)
            attributes(trd) <- attributes(y)
            if(trim) trd <- window(trd, start = index(y)[n.periods])
            return(trd)
          }
)
