#' Present value for geometric gradient series (Engineering Economics)
#'
#' Compute P given A1
#'
#' P is expressed as
#'
#' 	\deqn{P = \frac{A_1\left[1 - \left(1 + i\right)^{-n}\left(1 + f\right)^{n}\right]}{i - f}, \: where \: f \neq i}
#'
#' or
#'
#'   \deqn{P = A_1n\left(1 + i\right)^{-1}, \: where \: f = i}
#'
#' \describe{
#'	\item{\emph{P}}{"the present equivalent of the geometric gradient
#'     series"}
#'	\item{\emph{\eqn{A_1}}}{"the initial cash flow in that occurs at the end of
#'     period one"}
#'	\item{\emph{i}}{the "interest rate per period"}
#'	\item{\emph{f}}{the "average rate each period"}
#'   \item{\emph{n}}{the "number of interest periods"}
#' }
#'
#' 	Note: "f can be positive or negative"
#'
#'
#' @param A1 numeric vector that contains the initial annual value(s)
#' @param n numeric vector that contains the period value(s)
#' @param f numeric vector that contains the average interest rate value(s)
#'     as a percent per period
#' @param i numeric vector that contains the interest rate(s) as a percent
#'
#' @return PgivenA1 numeric vector that contains the present value(s) rounded
#'    to 2 decimal places
#'
#' @references
#' William G. Sullivan, Elin M. Wicks, and C. Patrick Koelling, \emph{Engineering Economy}, Fourteenth Edition, Upper Saddle River, New Jersey: Pearson/Prentice Hall, 2009, page 156-159.
#'
#' @encoding UTF-8
#'
#'
#' @examples
#' library(iemisc)
#' # Example 4-23 from the Reference text (page 158-159)
#' PgivenA1(1000, 25, 20, 4) # i is 25% and f is 20%
#'
#'
#' # Example 4-24 from the Reference text (page 159)
#' PgivenA1(1000, 25, -20, 4) # i is 25% and f is -20%
#'
#'
#' @export
PgivenA1 <- function (A1, i, f, n) {

i <- i / 100

f <- f / 100

if (f != i) {

PgivenA1 <- (A1 * (1 - (1 + i) ^ -n * (1 + f) ^ n)) / (i - f)
}

else {

PgivenA1 <- A1 * n * (1 + i) ^ -1
}

return(round(PgivenA1, digits = 2))
}
