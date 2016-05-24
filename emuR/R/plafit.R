##' Calculate the coefficients of a parabola
##' 
##' Fit a second ordered polynomial to a vector of values
##' 
##' The function fits a parabola (2nd order polynomial) following the method of
##' van Bergem, Speech Communication, 14, 1994, 143-162. The algorithm fixes
##' the parabola at the onset, midpoint, and offset of the vector i.e. such
##' htat the fitted parabola and original vector have the same values at these
##' points.
##' 
##' @param wav a vector or single column matrix of numeric values to which the
##' 2nd order polynomial is to be fitted.
##' @param fit if F, return the coefficients of the polynomial; if T, the
##' values of the polynomial are returned to the same length as the vector wav.
##' @param n in fitting the polynomial, linear time normalisation is first
##' applied to the input vector wav to 101 points. The polynomial is fitted
##' under the assumption that these points extend linearly in time between t =
##' -1 and t = 1 with t = 0 occurring at the temporal midpoint.
##' @return The function returns the coefficients of c0, c1, c2 in the parabola
##' y = c0 + c1t + c2t\eqn{\mbox{\textasciicircum}}{^}2 where t extends between
##' -1 and 1. The function can also be used to derive the values of the
##' parabola as a function of time from the coefficients.
##' @author Jonathan Harrington
##' @seealso \code{\link{dct}}
##' @keywords math
##' @examples
##' 
##' # fit a polynomial to a segment of fundamental frequency data
##' plafit(vowlax.fund[1,]$data)
##' 
##' # return the fitted values of the polynomial
##' plafit(vowlax.fund[1,]$data, fit=TRUE)
##' 
##' 
##' @export plafit
"plafit" <- function(wav, fit = FALSE, n = 101)
{
  if(!is.vector(wav) & !is.matrix(wav) )
    stop("input signal must be a vector or a one-columned matrix")
  if(is.matrix(wav) )
  {
    if(ncol(wav)!=1)
      stop("input signal must be a vector or a one-columned matrix")
  }
  if(is.vector(wav))
    nz <- names(wav)
  if(is.matrix(wav))
    nz <- dimnames(wav)[[1]]
  if(n %% 2 != 1)
    n <- n + 1
  N <- length(wav)
  a <- approx(wav, n=n)$y
  times <- seq(-1, 1, length=n)
  c0 <- a[times==0]
  c1 <- 0.5 * (a[n] - a[1])
  c2 <- 0.5 * (a[n] + a[1]) - c0
  if(fit)
  {
    y <- c0 + c1 * times + c2 * (times^2)
    result <- approx(y, n=N)$y
    names(result) <- nz
  }
  else
  {
    result <- c(c0, c1, c2)
    names(result) <- c("c0", "c1", "c2")
  }
  result
}
