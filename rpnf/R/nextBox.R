#' Determine the next box frontier for current quote(s) given a recent XO-status.
#' 
#' Note: offset should only be used for reversal calculation
#' @param quote  A single quote or a vector of quotes.
#' @param status A single character indicating the current XO-status.
#' @param boxsize A single numeric value, indicating the boxsize to be considered.
#' @param log TRUE, if logarithmic scales should be used.
#' @param offset A numeric value 
nextBox <- function(quote,status, boxsize=1, log=FALSE, offset=0) {
  if (!(is.numeric(quote))) {
    stop("Argument quote has to be numeric!")
  }
  if (!(is.character(status) & nchar(status)==1 )) {
    stop("Argument status has to be a character and of length 1!")
  }
  if (!(is.numeric(boxsize) & length(boxsize)==1)) {
    stop("Argument boxsize has to be numeric and of length 1!")
  }
  if(!(is.logical(log) & length(log))) {
    stop("Argument log has to be logical and of length 1!")
  }
  if (!(is.numeric(offset) & length(offset)==1)) {
    stop("Argument offset has to be numeric and of length 1!")
  }
  
  if (status == "X")
    box2upper(quote2box(quote=quote,boxsize=boxsize,log=log)+offset,boxsize=boxsize,log=log)
  else 
    box2lower(quote2box(quote=quote,boxsize=boxsize,log=log)+offset,boxsize=boxsize,log=log)
}
