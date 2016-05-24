coredata <- function(x, ...)
  UseMethod("coredata")

coredata.default <- function(x, ...) x

coredata.zoo <- function(x, ...)
{
  attr(x, "class") <- attr(x, "oclass")
  attr(x, "index") <- attr(x, "oclass") <- attr(x, "frequency") <- NULL
  return(x)
}

## # experimental coredata.zoo to take advantage of new C code contributed from xts
## .coredata.zoo <- function(x, ...) {
##   if(length(x) == 0)
##     return(vector(storage.mode(x)))
##   .Call("zoo_coredata", x, TRUE, PACKAGE = "zoo")  # second arg is to copy most attr, for compat with xts
## }

coredata.ts <- function(x, ...)
{
  x <- unclass(x)
  attr(x, "tsp") <- NULL
  return(x)
}

coredata.irts <- function(x, ...)
{
  return(x$value)
}

coredata.its <- function(x, ...)
{
  return(x@.Data)
}


"coredata<-" <- function(x, value)
{
  UseMethod("coredata<-")
}

"coredata<-.zoo" <- function(x, value)
{
  stopifnot(length(x) == length(value))
  if(!(is.vector(value) || is.factor(value) || is.matrix(value) || is.data.frame(value)))
    stop(paste(dQuote("value"), ": attempt to assign invalid coredata to zoo object"))
  if(is.matrix(value) || is.data.frame(value)) value <- as.matrix(value)
    
  x[] <- value  
  attr(x, "oclass") <- attr(value, "class")
  return(x)
}

"coredata<-.ts" <- function(x, value)
{
  stopifnot(length(x) == length(value))
  dim(value) <- dim(x)
  x[] <- value
  return(x)
}

"coredata<-.irts" <- function(x, value)
{
  stopifnot(length(x$value) == length(value))
  dim(value) <- dim(x$value)
  x$value[] <- value
  return(x)
}

"coredata<-.its" <- function(x, value)
{
  stopifnot(length(x@.Data) == length(value))
  dim(value) <- dim(x@.Data)
  x@.Data[] <- as.matrix(value)
  return(x)
}
