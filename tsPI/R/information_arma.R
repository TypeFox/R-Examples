#' Large Sample Approximation of Information Matrix for ARMA process
#'
#' Fortran implementation of \code{InformationMatrixARMA} function of
#' \code{FitARMA} package, except that the function uses the same
#' ARMA model definition as \code{arima}, where both the
#' AR and MA parts of the model are on the right side of the equation, i.e.
#' MA coefficients differ in sign compared to \code{InformationMatrixARMA}.
#'
#' @export
#' @param phi Autoregressive coefficients.
#' @param theta Moving average coefficients.
#' @return Large sample approximation of information matrix for ARMA process.
#' @references
#' \enumerate{
#' \item{Box, G. and Jenkins, G. (1970). Time Series Analysis: Forecasting and Control. San Francisco: Holden-Day.}
#'\item{McLeod, A. I. and Zhang, Y., (2007). Faster ARMA maximum likelihood estimation Computational Statistics & Data
#' Analysis 52(4) URL http://dx.doi.org/10.1016/j.csda.2007.07.020}
#' }
information_arma <- function(phi = NULL, theta = NULL){
  imat <- .Fortran(fapproxinfmat, length(phi),
    length(theta), as.double(c(phi, theta)), imat = diag(length(phi) + length(theta)))$imat
  imat[lower.tri(imat)] <- imat[upper.tri(imat)]
  imat
}
