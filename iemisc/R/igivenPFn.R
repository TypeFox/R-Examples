#' Interest rate given Future value, Number of periods, and Present value (Engineering Economics)
#'
#' Compute i given F, n, and P
#'
#' i is expressed as
#'
#' 	\deqn{i = \sqrt[n]{\frac{F}{P}} - 1}
#'
#' \describe{
#'	\item{\emph{i}}{the "effective interest rate per interest period"}
#'	\item{\emph{F}}{the "future equivalent"}
#'	\item{\emph{P}}{the "present equivalent"}
#'	\item{\emph{n}}{the "number of interest periods}
#' }
#'
#'
#' @param F numeric vector that contains the future value(s)
#' @param n numeric vector that contains the period value(s)
#' @param P numeric vector that contains the present value(s)
#'
#' @return i numeric vector that contains the effective interest rate as a
#' percent rounded to 2 decimal places
#'
#' @references
#' William G. Sullivan, Elin M. Wicks, and C. Patrick Koelling, \emph{Engineering Economy}, Fourteenth Edition, Upper Saddle River, New Jersey: Pearson/Prentice Hall, 2009, page 128-129, 142.
#'
#' @encoding UTF-8
#'
#'
#' @examples
#' library(iemisc)
#' # Example for equation 4-6 from the Reference text (page 128)
#' igivenPFn(P = 500, F = 1000, n = 10)
#'
#'
#'
#' @importFrom pracma nthroot
#'
#'
#'
#' @export
igivenPFn <- function (P, F, n) {

igivenPFn <- (nthroot(F / P, n) - 1) * 100 # interest rate as percent

return(round(igivenPFn, digits = 2))
}
