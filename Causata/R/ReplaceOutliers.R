ReplaceOutliers <- function(this, ...){
  # generic function to replace outliers
  UseMethod("ReplaceOutliers")
}


ReplaceOutliers.numeric <- function(this, lowerLimit=NULL, upperLimit=NULL, ...) {
  # Replaces outliers from a continuous variable x.
  # Values less than (<) the lowerLimit are replaced
  # with the lowerLimit, and values greater than (>) upperLimit are replaced with the 
  # upperLimit.
  # Missing values are ignored
  #
  # Author: Justin Hemann
  
  # check if lower limit is applied
  if (!is.null(lowerLimit)) {
    # a limit was provided, apply it
    idx <- this < lowerLimit
    if (any(is.na(idx))){
      idx[is.na(idx)] <- FALSE
    }
    this[idx] <- lowerLimit
  }
  
  # check if upper limit is applied
  if (!is.null(upperLimit)) {
    # a limit was provided, apply it
    idx <- this > upperLimit
    if (any(is.na(idx))){
      idx[is.na(idx)] <- FALSE
    }
    this[idx] <- upperLimit
  }
  return(this)
}