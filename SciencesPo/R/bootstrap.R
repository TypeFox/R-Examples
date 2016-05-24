#' @encoding UTF-8
#' @title Method for Bootstrapping
#' @description This method is intended to be provides statistical models that support  bootstrapping.
#' @param x is a vector or a fitted model object whose parameters will be used to produce bootstrapped statistics. Model objects are from the class \dQuote{glm} or \dQuote{lm}.
#' @param  \dots further arguments passed to or used by other methods.
#' @return A list with the \dQuote{alpha} and \dQuote{beta} slots set. Note that \dQuote{alpha} corresponds to ancillary parameters and \dQuote{beta} corresponds to systematic components of the model.
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#' @export
bootstrap <- function (x, ...)
  UseMethod("bootstrap")

#' @title Bootstrap
#'
#' @description This function is used to estimating standard errors when the distribution is not know.
#'
#' @param nboots The number of bootstraps.
#' @param FUN the name of the statistic to bootstrap, ie., 'mean', 'var', 'cov', etc as a string.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#'
#' @examples
#' x = runif(10, 0, 1)
#' bootstrap(x,FUN=mean)
#'
#' @rdname bootstrap
#' @export
bootstrap.default<-function(x, nboots=100, FUN,  ...){
	n=length(x)
	lings <-replicate(nboots, match.fun(FUN)(sample(x, n, replace=TRUE)))
	list(se = sd(lings),
       lings = lings)
}
NULL


#' @title Bootstrap Parameters of a Statistical Model
#'
#' @description This method is used to bootstrapping statistical models, typically of class \dQuote{lm} or \dQuote{glm}.
#'
#' @rdname bootstrap
#' @export
bootstrap.model <- function (x, ...)
  list(
    alpha = NULL,
    beta = coef(x)
  )
