##
##
## Generic for returnconvert
##
setGeneric(name = "returnconvert", def = function(y, convdir = c("cont2disc", "disc2cont"), percentage = TRUE){standardGeneric("returnconvert")})
##
## Methods for returnconvert
## (in alphabetical order)
##
## for class data.frame
setMethod(f = "returnconvert",
          signature = c(y = "data.frame"),
          definition = function(y, convdir = c("cont2disc", "disc2cont"), percentage = TRUE){
            if(!all(apply(y, 2, is.numeric))){
              stop("Data frame has non-numeric columns.\n")
            }
            r <- as.data.frame(apply(y, 2, returnconvert, convdir = convdir, percentage = percentage), ncol = ncol(y))
            colnames(r) <- colnames(y)
            return(r)
          }
)
## for class matrix
setMethod(f = "returnconvert",
          signature = c(y = "matrix"),
          definition = function(y, convdir = c("cont2disc", "disc2cont"), percentage = TRUE){
            r <- matrix(apply(y, 2, returnconvert, convdir = convdir, percentage = percentage), ncol = ncol(y))
            colnames(r) <- colnames(y)
            return(r)
          }
)
## for class mts
setMethod(f = "returnconvert",
          signature = c(y = "mts"),
          definition = function(y, convdir = c("cont2disc", "disc2cont"), percentage = TRUE){
            r <- matrix(apply(y, 2, returnconvert, convdir = convdir, percentage = percentage), ncol = ncol(y))
            attributes(r) <- attributes(y)
            return(r)
          }
)
## for class numeric
setMethod(f = "returnconvert",
          signature = c(y = "numeric"),
          definition = function(y, convdir = c("cont2disc", "disc2cont"), percentage = TRUE){
            convdir <- match.arg(convdir)
            r <- y
            if(percentage) r <- r / 100
            if(convdir == "cont2disc"){
              r <- exp(r) - 1.0
            }
            if(convdir == "disc2cont"){
              r <- log(1 + r)
            }
            if(percentage) r <- r * 100
            return(r)
          }
)
## for class timeSeries
setMethod(f = "returnconvert",
          signature = c(y = "timeSeries"),
          definition = function(y, convdir = c("cont2disc", "disc2cont"), percentage = TRUE){
            r <- apply(y, 2, returnconvert, convdir = convdir, percentage = percentage)
            r <- timeSeries(r, charvec = time(y))
            return(r)
          }
)
## for class ts
setMethod(f = "returnconvert",
          signature = c(y = "ts"),
          definition = function(y, convdir = c("cont2disc", "disc2cont"), percentage = TRUE){
            r <- returnconvert(c(y), convdir = convdir, percentage = percentage)
            attributes(r) <- attributes(y)
            return(r)
          }
)
## for class xts
setMethod(f = "returnconvert",
          signature = c(y = "xts"),
          definition = function(y, convdir = c("cont2disc", "disc2cont"), percentage = TRUE){
            yc <- as.matrix(coredata(y))
            r <- returnconvert(yc, convdir = convdir, percentage = percentage)
            attributes(r) <- attributes(y)
            return(r)
          }
)
## for class zoo
setMethod(f = "returnconvert",
          signature = c(y = "zoo"),
          definition = function(y, convdir = c("cont2disc", "disc2cont"), percentage = TRUE){
            yc <- as.matrix(coredata(y))
            r <- returnconvert(yc, convdir = convdir, percentage = percentage)
            attributes(r) <- attributes(y)
            return(r)
          }
)
