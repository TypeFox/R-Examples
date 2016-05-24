#' Determine an appropriate boxsize, if you want to use logarithmic scale.
#'
#' This function returns an appropriate boxsize if you want to do your point and figure analysis with an logarithmic scale.
#' 
#' @param percent a numeric value defining the percent 
#' @return a numeric value which is equivalent to the percental change given on a logarithmic scale
#' @export
#' @examples
#' # apply it with pnfprocessor
#' library(rpnf) # Load rpnf library
#' data(DOW) # Load some example data
#' 
#' # return appropriate value for 1% boxsize
#' getLogBoxsize(percent=1)
#' 
#' pnfprocessor(
#'  high=DOW$High,
#'  low=DOW$Low,
#'  date=DOW$Date,
#'  boxsize=getLogBoxsize(percent=1),
#'  log=TRUE)
getLogBoxsize <- function(percent) {
  if (!is.numeric(percent)) {
    stop("Argument percent has to be numeric!")
  }
  if (length(percent)>1) {
    warning("Argument percent for function getLogBoxsize() is a vector, and will return a vector!")
  }
  if (min(percent)<=0) {
    if (min(percent)<=-100) {
      warning("Argument percent contains values less equal than -100, this will introduce NaN's and strange results!")
    } else {
      warning("Argument percent contains values less equal than zero, this will introduce strange results!")
    }
  }
  return (log(100+percent)-log(100))
}