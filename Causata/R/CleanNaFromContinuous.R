
#
# define generic functions
#
CleanNaFromContinuous <- function(x, ...) {
  # generic function to replace missing values in a numeric variable
  UseMethod("CleanNaFromContinuous", x)
}


CleanNaFromContinuousFunction <- function(x, method="median", replacement.value=NULL, return.replacement=FALSE){
  label.median = "median"
  label.mean   = "mean"
  
  # replaces missing values in a vector
  if (!is.null(replacement.value)){ 
    # the replacement value was specified, so use it 
  } else if (method == label.median){
    replacement.value <- median( na.omit(x) )
  } else if (method == label.mean){
    replacement.value <- mean(   na.omit(x) )
  } else {
    stop("Invalid value of method provided, halting: ", method)
  }
  
  # replace missing values with the computed value from above
  x[is.na(x)] <- replacement.value
  
  # return the result
  if (return.replacement){
    # flag is true, so return a list with the vector and the replacement value
    return(list(x=x, replacement.value=replacement.value))
  } else {
    # return only the vector
    return(x)
  }
}

#
# class-specific functions
#

CleanNaFromContinuous.numeric <- function(x, method="median", replacement.value=NULL, return.replacement=FALSE, ...){
  return(CleanNaFromContinuousFunction(x, method, replacement.value, return.replacement))
}

CleanNaFromContinuous.POSIXct <- function(x, method="median", replacement.value=NULL, return.replacement=FALSE, ...){
  return(CleanNaFromContinuousFunction(x, method, replacement.value, return.replacement))
}