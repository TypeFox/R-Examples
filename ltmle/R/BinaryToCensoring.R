#' BinaryToCensoring
#' 
#' Helper function for creating censoring columns as factors.
#' 
#' Exactly one of \code{is.censored} and \code{is.uncensored} must be specified
#' as a \emph{named} argument.  All elements of the input vector must be 0, 1,
#' or \code{NA}
#' 
#' @param is.censored binary vector: 0=uncensored, 1=censored
#' @param is.uncensored binary vector: 0=censored, 1=uncensored
#' @return an object of class "\code{factor}" with levels "censored" and
#' "uncensored"
#' @author Joshua Schwab \email{joshuaschwab@@yahoo.com}
#' @seealso \code{\link{factor}}
#' @examples
#' 
#'  BinaryToCensoring(is.censored=c(0, 1, 1, 0, NA))
#'  BinaryToCensoring(is.uncensored=c(1, 0, 0, 1, NA))   #the same
#'  
#'  \dontrun{
#'  BinaryToCensoring(c(0, 1))   #error because the input must be named
#'  }
#' 
#' @export BinaryToCensoring
BinaryToCensoring <- function(is.censored, is.uncensored) {
  if (! xor(missing(is.censored), missing(is.uncensored))) stop("exactly one of is.censored and is.uncensored must be passed")
  calling.name <- names(sys.call(0))[2]
  if (length(calling.name) == 0 || ! calling.name %in% c("is.censored", "is.uncensored")) stop("the argument to BinaryToCensoring must be completely named - see ?BinaryToCensoring")
  if (missing(is.uncensored)) {
    is.uncensored <- ! is.censored
  }
  if (! all(is.uncensored %in% c(0, 1, NA))) stop("the argument to BinaryToCensoring should be binary (0, 1, or NA) or logical")
  y <- character(length(is.uncensored))
  y[is.uncensored == 0] <- "censored"
  y[is.uncensored == 1] <- "uncensored"
  y[is.na(is.uncensored)] <- NA
  return(factor(y))
}
