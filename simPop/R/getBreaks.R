#' Compute break points for categorizing (semi-)continuous variables
#'
#' Compute break points for categorizing continuous or semi-continuous
#' variables using (weighted) quantiles.  This is a utility function that is
#' useful for writing custom wrapper functions such as \code{\link{simEUSILC}}.
#'
#' If \code{equidist} is \code{TRUE}, the behavior is as follows.  If
#' \code{zeros} is \code{TRUE} as well, the 0\%, 10\%, \dots{}, 90\% quantiles
#' of the negative values and the 10\%, 20\%, \dots{}, 100\% of the positive
#' values are computed.  These quantiles are then used as break points together
#' with 0.  If \code{zeros} is not \code{TRUE}, on the other hand, the 0\%,
#' 10\%, \dots{}, 100\% quantiles of all values are used.
#'
#' If \code{equidist} is not \code{TRUE}, the behavior is as follows.  If
#' \code{zeros} is not \code{TRUE}, the 1\%, 5\%, 10\%, 20\%, 40\%, 60\%, 80\%,
#' 90\%, 95\% and 99\% quantiles of all values are used for the inner part of
#' the data (instead of the equidistant 10\%, \dots{}, 90\% quantiles).  If
#' \code{zeros} is \code{TRUE}, these quantiles are only used for the positive
#' values while the quantiles of the negative values remain equidistant.
#'
#' Note that duplicated values among the quantiles are discarded and that the
#' minimum and maximum are replaced with \code{lower} and \code{upper},
#' respectively, if these are specified.
#'
#' The (weighted) quantiles are computed with the function
#' \code{\link{quantileWt}}.
#'
#' @name getBreaks
#' @param x a numeric vector to be categorized.
#' @param weights an optional numeric vector containing sample weights.
#' @param zeros a logical indicating whether \code{x} is semi-continuous, i.e.,
#' contains a considerable amount of zeros. See \dQuote{Details} on how this
#' affects the behavior of the function.
#' @param lower,upper optional numeric values specifying lower and upper bounds
#' other than minimum and maximum of \code{x}, respectively.
#' @param equidist a logical indicating whether the (positive) break points
#' should be equidistant or whether there should be refinements in the lower
#' and upper tail (see \dQuote{Details}).
#' @param probs a numeric vector of probabilities with values in \eqn{[0, 1]}
#' giving quantiles to be used as (positive) break points.  If supplied, this
#' is preferred over \code{equidist}.
#' @param strata an optional vector specifying a strata variable (e.g household ids).
#' if specified, the mean of \code{x} (and also of \code{weights} if specified) is
#' computed within each strata before calculating the breaks.
#' @return A numeric vector of break points.
#' @author Andreas Alfons and Bernhard Meindl
#' @export
#' @seealso \code{\link{getCat}}, \code{\link{quantileWt}}
#' @keywords manip
#' @examples
#'
#' data(eusilcS)
#'
#' # semi-continuous variable, positive break points equidistant
#' getBreaks(eusilcS$netIncome, weights=eusilcS$rb050)
#'
#' # semi-continuous variable, positive break points not equidistant
#' getBreaks(eusilcS$netIncome, weights=eusilcS$rb050,
#'     equidist = FALSE)
#'
getBreaks <- function(x, weights = NULL, zeros = TRUE, lower = NULL,
  upper = NULL, equidist = TRUE, probs = NULL, strata=NULL) {

  . <- NULL

  if ( !is.null(strata) ) {
    if ( length(strata) != length(x) ) {
      stop("'strata' must have the same length as 'x'")
    }
    if ( is.null(weights) ) {
      dat <- data.table(x=x, strata=strata)
    } else {
      dat <- data.table(x=x, strata=strata, weights=weights)
    }
    setkey(dat, strata)
    if ( is.null(weights) ) {
      dat <- dat[,.(x=mean(x, na.rm=TRUE)), by=key(dat)]
    } else {
      dat <- dat[,.(x=mean(x, na.rm=TRUE), weights=mean(weights,na.rm=TRUE)), by=key(dat)]
      weights <- dat$weights
    }
    x <- dat$x
  }

  # initializations
  if(!is.numeric(x)) stop("'x' must be a numeric vector")
  if(!is.null(weights)) {
    if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
    else if(length(weights) != length(x)) {
      stop("'weights' must have the same length as 'x'")
    }
  }
  zeros <- isTRUE(zeros)
  if(!is.null(probs)) {
    if(!is.numeric(probs) || all(is.na(probs)) ||
       isTRUE(any(probs < 0 | probs > 1))) {
      stop("'probs' must be a numeric vector with values in [0,1]")
    }
  }
  if(zeros) {
    pos <- which(x > 0)
    if(length(pos)) {
      if(is.null(probs)) {
        if(isTRUE(equidist)) probs <- seq(0.1, 1, by=0.1)
        else probs <- c(0.01, 0.05, 0.1, seq(0.2, 0.8, by=0.2), 0.9, 0.95, 0.99, 1)
      } else probs <- c(probs, 1)
      qpos <- quantileWt(x[pos], weights[pos], probs)
    } else qpos <- NULL
    neg <- which(x < 0)
    if(length(neg)) {
      pneg <- seq(0, 0.9, by=0.1)
      qneg <- quantileWt(x[neg], weights[neg], pneg)
    } else qneg <- NULL
    breaks <- c(qneg, 0, qpos)
  } else {
    if(is.null(probs)) {
      if(isTRUE(equidist)) probs <- seq(0, 1, by=0.1)
      else probs <- c(0, 0.01, 0.05, 0.1, seq(0.2, 0.8, by=0.2), 0.9, 0.95, 0.99, 1)
    } else probs <- c(0, probs, 1)
    breaks <- quantileWt(x, weights, probs)
  }
  breaks <- unique(breaks)  # remove duplicated values
  if(!is.null(lower)) {
    if(!is.numeric(lower) || length(lower) == 0) {
      stop("'lower' is not numeric or has length 0")
    } else if(length(lower) > 1) lower <- lower[1]
    if(isTRUE(lower > breaks[1])) {
      warning("'lower' is larger than the smallest ",
              "breakpoint and therefore disregarded")
    } else breaks[1] <- lower
  }
  if(!is.null(upper)) {
    if(!is.numeric(upper) || length(upper) == 0) {
      stop("'upper' is not numeric or has length 0")
    } else if(length(upper) > 1) upper <- upper[1]
    nb <- length(breaks)
    if(isTRUE(upper < breaks[nb])) {
      warning("'upper' is smaller than the largest ",
              "breakpoint and therefore disregarded")
    } else breaks[nb] <- upper
  }
  breaks
}
