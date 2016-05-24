#' Relative error
#'
#' This function computes the relative error.
#'
#' Relative error is expressed as
#'
#' \deqn{\varepsilon_{t} = \frac{{true \: value} - {approximation}}{true \: value} \cdot 100}
#'
#' \describe{
#'	\item{\emph{\eqn{\varepsilon_t}}}{the "true percent relative error"}
#'	\item{\emph{true value}}{the true value}
#'	\item{\emph{approximation}}{the approximate value}
#' }
#'
#'
#' @param xt numeric vector that contains the true value(s)
#' @param xa numeric vector that contains the approximate value(s)
#'
#' @return relative error, as a percent (\%), as a numeric \code{\link{vector}}.
#'
#'
#' @references
#' Steven C. Chapra, \emph{Applied Numerical Methods with MATLAB for Engineers and Scientists}, Second Edition, Boston, Massachusetts: McGraw-Hill, 2008, page 82-83.
#'
#'
#' @encoding UTF-8
#'
#' @seealso \code{\link{sgm}} for geometric mean, \code{\link{shm}} for harmonic mean, \code{\link{cv}} for
#'  coefficient of variation (CV), \code{\link{rms}} for root-mean-square (RMS), \code{\link{approxerror}}
#'  for approximate error, and \code{\link{ranges}} for sample range.
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' @examples
#' library(iemisc)
#' # Example 4.1 from the Reference text (page 83)
#' relerror(1.648721, 1.5) # answer as a percent (\%)
#'
#'
#'
#' @export
relerror <- function (xt, xa) {

abs(((xt - xa) / xt)) * 100

}




#' Approximate error
#'
#' This function computes the "approximate estimate of the error" ("percent relative error").
#'
#' Approximate error is expressed as
#'
#' \deqn{\varepsilon_{a} = \frac{{present \: approximation} - {previous \: approximation}}{present \: approximation} \cdot 100}
#'
#' \describe{
#'	\item{\emph{\eqn{\varepsilon_a}}}{the "approximate estimate of the error"}
#'	\item{\emph{present approximation}}{the "present approximation"}
#'	\item{\emph{previous approximation}}{the "previous approximation"}
#' }
#'
#'
#' @param pres numeric vector that contains the "present approximation"
#'   value(s)
#' @param prev numeric vector that contains the "previous approximation"
#'   value(s)
#'
#' @return approximate error, as a percent (\%), as a numeric \code{\link{vector}}.
#'
#'
#' @references
#' Steven C. Chapra, \emph{Applied Numerical Methods with MATLAB for Engineers and Scientists}, Second Edition, Boston, Massachusetts: McGraw-Hill, 2008, page 82-84.
#'
#'
#' @encoding UTF-8
#'
#' @seealso \code{\link{sgm}} for geometric mean, \code{\link{shm}} for harmonic mean, \code{\link{cv}} for
#'  coefficient of variation (CV), \code{\link{rms}} for root-mean-square (RMS), \code{\link{relerror}}
#'  for relative error, and \code{\link{ranges}} for sample range.
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' @examples
#' library(iemisc)
#' # Example 4.1 from the Reference text (page 84)
#' approxerror(1.5, 1) # answer as a percent (\%)
#'
#'
#'
#' @export
approxerror <- function (pres, prev) {

abs(((pres - prev) / pres)) * 100

}
