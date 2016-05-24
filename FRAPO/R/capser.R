##
##
## Generic for capser trend 
##
setGeneric(name = "capser", def = function(y, min, max){standardGeneric("capser")})
##
## Methods for capser
## (in alphabetical order)
##
## for class data.frame
setMethod(f = "capser",
          signature = c(y = "data.frame"),
          definition = function(y, min, max){
            if(!all(apply(y, 2, is.numeric))){
              stop("Data frame has non-numeric columns.\n")
            }
            cs <- data.frame(apply(y, 2, capser, min = min, max = max))
            return(cs)
          }
)
## for class matrix
setMethod(f = "capser",
          signature = c(y = "matrix"),
          definition = function(y, min, max){
            cs <- apply(y, 2, capser, min = min, max = max)
            return(cs)
          }
)
## for class mts
setMethod(f = "capser",
          signature = c(y = "mts"),
          definition = function(y, min, max){
            cs <- apply(y, 2, capser, min = min, max = max)
            attributes(cs) <- attributes(y)
            return(cs)
          }
)
## for class numeric
setMethod(f = "capser",
          signature = c(y = "numeric"),
          definition = function(y, min, max){
            min <- as.numeric(min)[1]
            max <- as.numeric(max)[1]
            if(min >= max){
              stop("\nMinimum value is greater or equal than for maximum.\n")
            }
            cs <- y
            cs[y < min] <- min
            cs[y > max] <- max
            return(cs)
          }
)
## for class timeSeries
setMethod(f = "capser",
          signature = c(y = "timeSeries"),
          definition = function(y, min, max){
            cs <- apply(y, 2, capser, min = min, max = max)
            return(cs)
          }
)
## for class ts
setMethod(f = "capser",
          signature = c(y = "ts"),
          definition = function(y, min, max){
            cs <- capser(c(y), min = min, max = max)
            attributes(cs) <- attributes(y)
            return(cs)
          }
)
## for class xts
setMethod(f = "capser",
          signature = c(y = "xts"),
          definition = function(y, min, max){
            yc <- as.matrix(coredata(y))
            cs <- apply(yc, 2, capser, min = min, max = max)
            attributes(cs) <- attributes(y)
            return(cs)
          }
)
## for class zoo
setMethod(f = "capser",
          signature = c(y = "zoo"),
          definition = function(y, min, max){
            yc <- as.matrix(coredata(y))
            cs <- apply(yc, 2, capser, min = min, max = max)            
            attributes(cs) <- attributes(y)
            return(cs)
          }
)
