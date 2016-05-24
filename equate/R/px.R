#' Percentile Ranks and Cumulativie Frequencies
#' 
#' These functions compute percentile ranks and cumulative frequency
#' distributions for frequency tables.
#' 
#' These functions compute percentile ranks and cumulative frequencies for a
#' univariate distribution, and percentile ranks from one univariate
#' distribution (\code{x}) corresponding to score values in another (\code{y}).
#' 
#' @param x either a vector of counts, or an object of class
#' \dQuote{\code{freqtab}} from which counts will be taken.
#' @param y an object of class \dQuote{\code{freqtab}} when \code{x} is as
#' well, otherwise, a vector or \code{data.frame} of counts. See below for
#' details.
#' @param ys vector specifying the \code{y} score scale, when it is not
#' contained in the first column of \code{y}. If \code{y} can be converted to a
#' \code{data.frame}, it is assumed to be univariate with the first column
#' containing the score scale and the second containing the counts.
#' @param margin,ymargin integers specifying the margins for which frequencies
#' or percentile ranks will be returned. \code{margin} applies to \code{x} and
#' \code{ymargin} to \code{y}.
#' @param \dots further arguments passed to or from other methods.
#' @return A vector is returned containing either percentile ranks or
#' cumulative frequencies with length equal to \code{length(x)}.
#' @author Anthony Albano \email{tony.d.albano@@gmail.com}
#' @seealso \code{\link{freqtab}}
#' @keywords univar
#' @examples
#' 
#' x <- as.freqtab(ACTmath[, 1:2], drop = TRUE)
#' y <- as.freqtab(ACTmath[, c(1, 3)], drop = TRUE)
#' 
#' # Percentile ranks for the x scale
#' round(px(x), 3)
#' 
#' # Percentile ranks in y for x each score
#' round(px(x, y = y), 3)
#' 
#' # Cumulative frequency distribution for x
#' round(fx(x), 3)
#' 
#' @export px
px <- function(x, ...) UseMethod("px")

# @describeIn px Default precentile rank method for a vector of frequencies.
#' @rdname px
#' @export
px.default <- function(x, y, ys, ...) {
  
  if (missing(y)) {
    x <- as.numeric(x/sum(x))
    p <- .5 * x[1]
    for (i in 2:length(x))
      p[i] <- sum(x[1:(i - 1)]) + .5 * x[i]
  }
  else {
    y <- as.data.frame(y)
    if (ncol(y) == 2) {
      ys <- y[, 1]
      y <- y[, 2]
    }
    xs <- floor(x + .5)
    yn <- sum(y)
    f <- sapply(xs, function(xi)
      sum(y[ys <= xi])/yn)
    flow <- sapply(xs, function(xi)
      sum(y[ys <= xi - 1])/yn)
    p <- flow + (x - (xs - .5)) * (f - flow)
  }
  return(p)
}

# @describeIn px Percentile rank method for \dQuote{\code{\link{freqtab}}}
# objects.
#' @rdname px
#' @export
px.freqtab <- function(x, margin = 1,
  y, ymargin = 1, ...) {
  
  if (!margin %in% seq(margins(x)))
    stop("'margin' not found in 'x'")
  if (missing(y))
    p <- px.default(margin(x, margin))
  else {
    if (!ymargin %in% seq(margins(y)))
      stop("'ymargin' not found in 'y'")
    p <- px.default(scales(margin(x, margin)),
      as.data.frame(margin(y, ymargin)))
  }
  
  return(p)
}

#' @rdname px
#' @export
fx <- function(x, ...) UseMethod("fx")

# @describeIn px Default cumulative frequency distribution method for
# a vector of frequencies.
#' @rdname px
#' @export
fx.default <- function(x, ...) {
  
  return(as.numeric(cumsum(x/sum(x))))
}

# @describeIn px Cumulative frequency distribution method for
# \dQuote{\code{\link{freqtab}}} objects.
#' @rdname px
#' @export
fx.freqtab <- function(x, margin = 1, ...) {
  
  if (!margin %in% seq(margins(x)))
    stop("'margin' not found in 'x'")
  
  return(fx.default(margin(x, margin)))
}
