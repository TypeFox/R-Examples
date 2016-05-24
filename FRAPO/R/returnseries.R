##
##
## Generic for returnseries
##
setGeneric(name = "returnseries", def = function(y, method = c("continuous", "discrete"), percentage = TRUE, trim = FALSE, compound = FALSE){standardGeneric("returnseries")}) 
##
## Methods for returnseries
## (in alphabetical order)
##
## for class data.frame
setMethod(f = "returnseries",
          signature = c(y = "data.frame"),
          definition = function(y, method = c("continuous", "discrete"), percentage = TRUE, trim = FALSE, compound = FALSE){
            if(!all(apply(y, 2, is.numeric))){
              stop("Data frame has non-numeric columns.\n")
            }
            ret <- as.data.frame(apply(y, 2, returnseries, method = method, percentage = percentage, trim = trim, compound = compound), ncol = ncol(y))
            colnames(ret) <- colnames(y)
            ifelse(trim, rownames(ret) <- rownames(y)[-1], rownames(ret) <- rownames(y)) 
            return(ret)
          }
)
## for class matrix
setMethod(f = "returnseries",
          signature = c(y = "matrix"),
          definition = function(y, method = c("continuous", "discrete"), percentage = TRUE, trim = FALSE, compound = FALSE){
            ret <- matrix(apply(y, 2, returnseries, method = method, percentage = percentage, trim = trim, compound = compound), ncol = ncol(y))
            colnames(ret) <- colnames(y)
            return(ret)
          }
)
## for class mts
setMethod(f = "returnseries",
          signature = c(y = "mts"),
          definition = function(y, method = c("continuous", "discrete"), percentage = TRUE, trim = FALSE, compound = FALSE){
            ret <- matrix(apply(y, 2, returnseries, method = method, percentage = percentage, trim = FALSE, compound = compound), ncol = ncol(y))
            attributes(ret) <- attributes(y)
            if(trim) ret <- window(ret, start = time(y)[2])
            return(ret)
          }
)
## for class numeric
setMethod(f = "returnseries",
          signature = c(y = "numeric"),
          definition = function(y, method = c("continuous", "discrete"), percentage = TRUE, trim = FALSE, compound = FALSE){
            y <- na.fail(y)
            method <- match.arg(method)
            if(method == "continuous"){
              ret <- c(NA, diff(log(y)))
              if(compound){
                ret[1] <- 0
                ret <- cumsum(ret)
              }
            }
            if(method == "discrete"){
              ret <- c(NA, diff(y) / y[-length(y)])
              if(compound){
                ret[1] <- 0
                ret <- cumprod(1 + ret) - 1
              }
            }
            if(percentage) ret <- ret * 100
            if(trim) ret <- ret[-1]
            return(ret)
          }
)
## for class timeSeries
setMethod(f = "returnseries",
          signature = c(y = "timeSeries"),
          definition = function(y, method = c("continuous", "discrete"), percentage = TRUE, trim = FALSE, compound = FALSE){
            ret <- apply(y, 2, returnseries, method = method, percentage = percentage, trim = FALSE, compound = compound)
            ret <- timeSeries(ret, charvec = time(y))
            if(trim) ret <- window(ret, start = time(y)[2], end = end(y)) 
            return(ret)
          }
)
## for class ts
setMethod(f = "returnseries",
          signature = c(y = "ts"),
          definition = function(y, method = c("continuous", "discrete"), percentage = TRUE, trim = FALSE, compound = FALSE){
            ret <- returnseries(c(y), method = method, percentage = percentage, trim = FALSE, compound = compound)
            attributes(ret) <- attributes(y)
            if(trim) ret <- window(ret, start = time(y)[2])
            return(ret)
          }
)
## for class xts
setMethod(f = "returnseries",
          signature = c(y = "xts"),
          definition = function(y, method = c("continuous", "discrete"), percentage = TRUE, trim = FALSE, compound = FALSE){
            yc <- as.matrix(coredata(y))
            ret <- returnseries(yc, method = method, percentage = percentage, trim = FALSE, compound = compound)
            attributes(ret) <- attributes(y)
            if(trim) ret <- window(ret, start = index(y)[2])
            return(ret)
          }
)
## for class zoo
setMethod(f = "returnseries",
          signature = c(y = "zoo"),
          definition = function(y, method = c("continuous", "discrete"), percentage = TRUE, trim = FALSE, compound = FALSE){
            yc <- as.matrix(coredata(y))
            ret <- returnseries(yc, method = method, percentage = percentage, trim = FALSE, compound = compound)
            attributes(ret) <- attributes(y)
            if(trim) ret <- window(ret, start = index(y)[2])
            return(ret)
          }
)
