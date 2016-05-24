#' Simple Interest Paid (Engineering Economics)
#'
#' Computes the total amount paid at the end of n periods using simple interest
#'
#' Simple Interest is expressed as
#'
#' \deqn{I = Pni}
#'
#' \deqn{S_n = P + I}
#'
#' or
#'
#' \deqn{S_n = P\left(1 + ni\right)}
#'
#' \describe{
#'	\item{\emph{P}}{the "principal amount (lent or borrowed)"}
#'	\item{\emph{\eqn{S_n}}}{the "total amount paid back"}
#'	\item{\emph{I}}{the "simple interest"}
#'	\item{\emph{i}}{the "interest rate per interest period"}
#'	\item{\emph{n}}{the "number of interest periods"}
#' }
#'
#'
#' @param P numeric vector that contains the present value(s)
#' @param n numeric vector that contains the period value(s)
#' @param i numeric vector that contains the interest rate(s) as whole number
#'    or decimal
#'
#' @return SimpIntPaid numeric vector that contains the total amount paid at
#'    the end of n periods rounded to 2 decimal places
#'
#' @references
#' \enumerate{
#' \item Chinyere Onwubiko, \emph{An Introduction to Engineering}, Mission, Kansas: Schroff Development Corporation, 1997, page 205-206.
#' \item William G. Sullivan, Elin M. Wicks, and C. Patrick Koelling, \emph{Engineering Economy}, Fourteenth Edition, Upper Saddle River, New Jersey: Pearson/Prentice Hall, 2009, page 116.
#' }
#'
#' @encoding UTF-8
#'
#'
#' @examples
#' library(iemisc)
#' # Example for equation 4-1 from the Sullivan Reference text (page 116)
#' SimpIntPaid(1000, 3, 10) # the interest rate is 10%
#'
#'
#'
#'
#' @export
SimpIntPaid <- function (P, n, i) {

i <- i / 100

I <- P * n * i # total interest using simple interest

return(round(P * (1 + n * i), digits = 2)) # total amount paid
}



#' Compound Interest Paid (Engineering Economics)
#'
#' Computes the total amount paid at the end of n periods using compound
#' interest
#'
#' Compound Interest is expressed as
#'
#' \deqn{S_n = P\left(1 + i\right)^n}
#'
#' \describe{
#'	\item{\emph{P}}{the "principal amount (lent or borrowed)"}
#'	\item{\emph{\eqn{S_n}}}{the "total amount paid back"}
#'	\item{\emph{i}}{the "interest rate per interest period"}
#'	\item{\emph{n}}{the "number of interest periods"}
#' }
#'
#'
#' @param P numeric vector that contains the present value(s)
#' @param n numeric vector that contains the period value(s)
#' @param i numeric vector that contains the interest rate(s) as a percent
#' @param frequency character vector that contains the frequency used to
#'    obtain the number of periods [annual (1), semiannual (2), quarter (4),
#'    bimonth (6), month (12), daily (365)]
#'
#' @return CompIntPaid numeric vector that contains the total amount paid at
#'    the end of n periods rounded to 2 decimal places
#'
#' @references
#' \enumerate{
#'  \item \emph{SFPE Handbook of Fire Protection Engineering}. 3rd Edition, DiNenno, P. J.; Drysdale, D.; Beyler, C. L.; Walton, W. D., Editor(s), 5/93-104 p., 2002. Chapter 7; Section 5; NFPA HFPE-02. See \url{http://fire.nist.gov/bfrlpubs//build02/art155.html}.
#'  \item William G. Sullivan, Elin M. Wicks, and C. Patrick Koelling, \emph{Engineering Economy}, Fourteenth Edition, Upper Saddle River, New Jersey: Pearson/Prentice Hall, 2009, page 120.
#' \item Chinyere Onwubiko, \emph{An Introduction to Engineering}, Mission, Kansas: Schroff Development Corporation, 1997, page 205-206.
#' }
#'
#'
#' @encoding UTF-8
#'
#'
#' @examples
#' library(iemisc)
#' # Compound Interest example from SFPE Reference text
#' CompIntPaid(100, 5, 10, frequency = "annual") # the interest rate is 10%
#'
#'
#'
#'
#' @export
CompIntPaid <- function (P, n, i, frequency = c("annual", "semiannual", "quarter", "bimonth", "month", "daily")) {

i <- i / 100

fr <- frequency

if (fr == "annual") {
fr <- 1
n <- n * fr

i <- i / fr

CompIntPaid <- P * (1 + i) ^ n # total amount paid

return(round(CompIntPaid, digits = 2))

} else if (fr == "semiannual") {

fr <- 2
n <- n * fr

i <- i / fr

CompIntPaid <- P * (1 + i) ^ n # total amount paid

return(round(CompIntPaid, digits = 2))

} else if (fr == "quarter") {

fr <- 4
n <- n * fr

i <- i / fr

CompIntPaid <- P * (1 + i) ^ n # total amount paid

return(round(CompIntPaid, digits = 2))

} else if (fr == "bimonth") {

fr <- 6
n <- n * fr

i <- i / fr

CompIntPaid <- P * (1 + i) ^ n # total amount paid

return(round(CompIntPaid, digits = 2))

} else if (fr == "month") {

fr <- 12
n <- n * fr

i <- i / fr

CompIntPaid <- P * (1 + i) ^ n # total amount paid

return(round(CompIntPaid, digits = 2))

} else if (fr == "daily") {

fr <- 365
n <- n * fr

i <- i / fr

CompIntPaid <- P * (1 + i) ^ n # total amount paid

return(round(CompIntPaid, digits = 2))

}
}
