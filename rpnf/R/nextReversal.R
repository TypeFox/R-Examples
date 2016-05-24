#' Determine the next reversal frontier for current quote(s) given a recent XO-status.
#' @param quote  A single quote or a vector of quotes.
#' @param status A single character indicating the current XO-status.
#' @param reversal number of boxes needed to make a reversal 
#' @param boxsize A single numeric value, indicating the boxsize to be considered.
#' @param log TRUE, if logarithmic scales should be used.
nextReversal <- function(quote,status, reversal=3L, boxsize=1, log=FALSE) {
  if (!(is.numeric(quote))) {
    stop("Argument quote has to be numeric!")
  }
  if (!(is.character(status) & nchar(status)==1)) {
    stop("Argument status has to be a character and of length 1!")
  }
  if (status == "X")
    return (nextBox(quote, "O", boxsize=boxsize,log=log, offset=-reversal+1))
  else
    return (nextBox(quote, "X", boxsize=boxsize,log=log, offset=reversal-1))
}
