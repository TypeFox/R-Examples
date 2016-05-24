#' @name confint
#' @author Sven E. Templer
#' @title Confidence intervals for numeric vectors
#' @description
#' Calculate confidence intervals for values of a numeric vector.
#' @param object A numeric vector.
#' @param parm Function for quantile calculation.
#' e.g. \code{qnorm}, \code{qt}
#' @param level Size of confidence (0 < size < 1).
#' @param ... Unused.
#' @param na.rm Logical, remove missing values for \code{sd} and \code{mean}.
#' @param ret.attr Logical, to include the mean value and function arguments
#' as attributes of the returned object.
#' @return 
#' Returns a numeric vector with the lower and upper range of the
#' confidence interval.
#' @examples
#' #
#' 
#' confint(1:3)
#' confint(1:3, ret.attr = FALSE)
#' 
#' #

#' @rdname confint
#' @export
confint.numeric <- function (object, parm = qnorm, level = .95, ..., 
                             na.rm = TRUE, ret.attr = TRUE) {
  if (any(level >= 1, level <= 0)) stop("level out of range")
  level <- 1 - (1 - level)/2
  n <- length(object)
  s <- sd(object, na.rm = na.rm)
  r <- parm(level) * s / sqrt(n)
  if (ret.attr) {
    m <- mean(object, na.rm = na.rm)
    mr <- c(m - r, m + r)
    attr(r, "mean") <- m
    attr(r, "range") <- mr
    attr(r, "level") <- level
    attr(r, "quantile") <- parm
  }
  return(r)
}