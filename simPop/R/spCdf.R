#' (Weighted empirical) cumulative distribution function
#' 
#' Compute a (weighted empirical) cumulative distribution function for survey
#' or population data.  For survey data, sample weights are taken into account.
#' 
#' Sample weights are taken into account by adjusting the step height.  To be
#' precise, the weighted step height for an observation is defined as its
#' weight divided by the sum of all weights\eqn{\ ( w_{i} / \sum_{j = 1}^{n}
#' w_{j} ).}{.}
#' 
#' If requested, the approximation is performed using the function
#' \code{\link[stats:approxfun]{approx}}.
#' 
#' @name spCdf
#' @param x a numeric vector.
#' @param weights an optional numeric vector containing sample weights.
#' @param approx a logical indicating whether an approximation of the
#' cumulative distribution function should be computed.
#' @param n a single integer value; if \code{approx} is \code{TRUE}, this
#' specifies the number of points at which the approximation takes place (see
#' \code{\link[stats:approxfun]{approx}}).
#' @return A list of class \code{"spCdf"} with the following components:
#' \item{x}{a numeric vector containing the \eqn{x}-coordinates.} \item{y}{a
#' numeric vector containing the \eqn{y}-coordinates.} \item{approx}{a logical
#' indicating whether the coordinates represent an approximation.}
#' @author Andreas Alfons and Stefan Kraft
#' @seealso \code{\link{spCdfplot}}, \code{\link[stats]{ecdf}},
#' \code{\link[stats:approxfun]{approx}}
#' @keywords dplot
#' @export
#' @examples
#' 
#' data(eusilcS)
#' cdfS <- spCdf(eusilcS$netIncome, weights = eusilcS$rb050)
#' plot(cdfS, type="s")
#' 
spCdf <- function(x, weights = NULL, approx = FALSE, n = 10000) {
  ## initializations
  # remove non-finites values
  ok <- is.finite(x)
  x <- x[ok]
  m <- length(x)
  if ( m == 0 ) {
    stop("'x' must contain finite values!\n")
  }
  # order observations
  ord <- sort.list(x, na.last=NA, method = "quick")
  x <- x[ord]
  # check whether CDF should be approximated
  approx <- isTRUE(approx)
  if ( approx ) {
    if ( !is.numeric(n) && length(n) != 1 && n <= 0 ) {
      stop("'n' must be a single positive integer")
    } else if(m < n) {
      approx <- FALSE
      warning("number of finite values in 'x' is smaller than 'n': no approximation")
    }
  }
  # define coordinates
  if ( is.null(weights) ) {
    y <- (1:m)/m
  } else {
    weights <- weights[ok][ord]
    cw <- cumsum(weights)
    y <- cw / cw[m]
  }
  res <- if(approx) approx(x, y, ties="ordered", n=n) else list(x=x, y=y)
  res$approx <- approx
  class(res) <- "spCdf"
  invisible(res)
}
