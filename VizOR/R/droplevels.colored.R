##' Drop unused levels from a colored factor and its color key
##'
##' This subclass method extends \code{droplevels.factor} to handle
##' updating of the special color key attribute of colored factors.
##' @title Drop unused levels from a colored factor
##' @param x a colored factor
##' @param ... further arguments passed to methods
##' @return a colored factor with no unused levels
##' @method droplevels colored
##' @S3method droplevels colored
##' @author David C. Norris
droplevels.colored <- function(x, ...){
  color.key <- key(x)
  x <- droplevels.factor(x)
  colored(x, color.key[levels(x)])
}
