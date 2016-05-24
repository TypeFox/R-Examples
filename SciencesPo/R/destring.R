#' @encoding UTF-8
#' @title Factors to numeric
#'
#' @description Converts factors to numeric like in Stata.
#'
#' @param x A factor whose levels will be converted.
#'
#' @seealso \code{\link{safe.chars}}.
#' @examples
#' mylevels <- c('Strongly Disagree', 'Disagree', 'Neither', 'Agree', 'Strongly Agree')
#'
#' myvar <- factor(sample(mylevels[1:5], 10, replace=TRUE))
#'
#' unclass(myvar) # testing order
#'
#' destring(myvar)
#'
#' @keywords Manipulation
#'
#' @export
`destring` <- function(x) {
  ## convert factor to strings
  if(is.character(x)) {
    as.numeric(x)
  } else if (is.factor(x)) {
    as.numeric(as.factor(x))
  } else if (is.numeric(x)) {
    invisible(x)
  } else {
    stop("Could not convert to numeric")
  }}
NULL
