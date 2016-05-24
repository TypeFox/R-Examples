#' Converts a single or a vector of quotes into integer boxnumbers for P&F-Analysis.
#' 
#' @param quote a single quote, or a vector of quotes
#' @param boxsize single numeric value, used as the boxsize
#' @param log TRUE, if logarithmic scales should be used
#' @return a single or a vector of integer boxnumbers
#' This function transforms a given quote into an unique integer box number
quote2box <- function(quote, boxsize=1, log=FALSE) {
  if (!is.numeric(quote)) {
    stop("Argument quote has to be numeric!")
  }
  if (!is.numeric(boxsize)) {
    stop("Argument boxsize has to be numeric!")
  }
  if (!is.logical(log)) {
    stop("Argument log has to be logical")
  }
  if (log & min(quote)<=0) {
    stop("Argument quotes must be greater than zero, if log=TRUE!")
  }
  if (length(boxsize)>1){
    stop("Argument boxsize is vector of length greater than 1. This is not supported yet!")
  }
  
  if (log==TRUE) {
    mylog <- function(x) {
      log(x)
    }
  } else {
    mylog <- function(x) {
      x
    }
  }
  result <- as.integer(floor(mylog(quote)/boxsize))
  result
}
