#' Annual value given Future value (Engineering Economics)
#'
#' Compute A given F
#'
#' A is expressed as
#'
#' 	\deqn{A = F\left[\frac{i}{\left(1 + i\right)^n - 1}\right]}
#'
#' \describe{
#'	\item{\emph{A}}{the "uniform series amount (occurs at the end of each
#'     interest period)"}
#'	\item{\emph{F}}{the "future equivalent"}
#'	\item{\emph{i}}{the "effective interest rate per interest period"}
#'	\item{\emph{n}}{the "number of interest periods"}
#' }
#'
#'
#' @param F numeric vector that contains the future value(s)
#' @param n numeric vector that contains the period value(s)
#' @param i numeric vector that contains the interest rate(s) as a percent
#' @param frequency character vector that contains the frequency used to
#'    obtain the number of periods [annual (1), semiannual (2), quarter (4),
#'    bimonth (6), month (12), daily (365)]
#'
#' @return AgivenF numeric vector that contains the annual value(s) rounded to
#'    2 decimal places
#' @return AF data.frame of both n (0 to n) and the resulting annual values
#'    rounded to 2 decimal places
#'
#' @references
#' William G. Sullivan, Elin M. Wicks, and C. Patrick Koelling, \emph{Engineering Economy}, Fourteenth Edition, Upper Saddle River, New Jersey: Pearson/Prentice Hall, 2009, page 135-136, 142, 164.
#'
#' @encoding UTF-8
#'
#'
#' @examples
#' library(iemisc)
#' # Example for equation 4-12 from the Reference text (page 135-136)
#' AgivenF(309*10^6, 60, 0.5, "annual")
#' # the interest rate is 0.5% per month and n is 60 months
#'
#' AF(309*10^6, 60, 0.5, "annual")
#' # the interest rate is 0.5% per month and n is 60 months
#'
#'
#' @import data.table
#'
#' @name AgivenF
NULL

#' @export
#' @rdname AgivenF
AgivenF <- function (F, n, i, frequency = c("annual", "semiannual", "quarter", "bimonth", "month", "daily")) {

i <- i / 100

fr <- frequency

if (fr == "annual") {
fr <- 1
n <- n * fr

i <- i / fr

AgivenF <- F * (i / (((1 + i) ^ n) - 1))

return(round(AgivenF, digits = 2))

} else if (fr == "semiannual") {

fr <- 2
n <- n * fr

i <- i / fr

AgivenF <- F * (i / (((1 + i) ^ n) - 1))

return(round(AgivenF, digits = 2))

} else if (fr == "quarter") {

fr <- 4
n <- n * fr

i <- i / fr

AgivenF <- F * (i / (((1 + i) ^ n) - 1))

return(round(AgivenF, digits = 2))

} else if (fr == "bimonth") {

fr <- 6
n <- n * fr

i <- i / fr

AgivenF <- F * (i / (((1 + i) ^ n) - 1))

return(round(AgivenF, digits = 2))

} else if (fr == "month") {

fr <- 12
n <- n * fr

i <- i / fr

AgivenF <- F * (i / (((1 + i) ^ n) - 1))

return(round(AgivenF, digits = 2))

} else if (fr == "daily") {

fr <- 365
n <- n * fr

i <- i / fr

AgivenF <- F * (i / (((1 + i) ^ n) - 1))

return(round(AgivenF, digits = 2))

}
}



#' @export
#' @rdname AgivenF
AF <- function (F, n, i, frequency = c("annual", "semiannual", "quarter", "bimonth", "month", "daily")) {

i <- i / 100

fr <- frequency

if (fr == "annual") {
fr <- 1
n <- n * fr

AF <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
AF[[y]] <- F * (i / (((1 + i) ^ seq(n)) - 1))
}

AF <- data.table(seq(n), unlist(AF))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

AF <- rbind(P0, AF)

setnames(AF, c("n (periods)", "Annual Worth ($US)"))
AF <- setDF(AF)

return(round(AF, digits = 2))

} else if (fr == "semiannual") {

fr <- 2
n <- n * fr

i <- i / fr

AF <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
AF[[y]] <- F * (i / (((1 + i) ^ seq(n)) - 1))
}

AF <- data.table(seq(n), unlist(AF))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

AF <- rbind(P0, AF)

setnames(AF, c("n (periods)", "Annual Worth ($US)"))
AF <- setDF(AF)

return(round(AF, digits = 2))

} else if (fr == "quarter") {

fr <- 4
n <- n * fr

i <- i / fr

AF <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
AF[[y]] <- F * (i / (((1 + i) ^ seq(n)) - 1))
}

AF <- data.table(seq(n), unlist(AF))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

AF <- rbind(P0, AF)

setnames(AF, c("n (periods)", "Annual Worth ($US)"))
AF <- setDF(AF)

return(round(AF, digits = 2))

} else if (fr == "bimonth") {

fr <- 6
n <- n * fr

i <- i / fr

AF <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
AF[[y]] <- F * (i / (((1 + i) ^ seq(n)) - 1))
}

AF <- data.table(seq(n), unlist(AF))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

AF <- rbind(P0, AF)

setnames(AF, c("n (periods)", "Annual Worth ($US)"))
AF <- setDF(AF)

return(round(AF, digits = 2))

} else if (fr == "month") {

fr <- 12
n <- n * fr

i <- i / fr

AF <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
AF[[y]] <- F * (i / (((1 + i) ^ seq(n)) - 1))
}

AF <- data.table(seq(n), unlist(AF))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

AF <- rbind(P0, AF)

setnames(AF, c("n (periods)", "Annual Worth ($US)"))
AF <- setDF(AF)

return(round(AF, digits = 2))

} else if (fr == "daily") {

fr <- 365
n <- n * fr

i <- i / fr

AF <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
AF[[y]] <- F * (i / (((1 + i) ^ seq(n)) - 1))
}

AF <- data.table(seq(n), unlist(AF))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

AF <- rbind(P0, AF)

setnames(AF, c("n (periods)", "Annual Worth ($US)"))
AF <- setDF(AF)

return(round(AF, digits = 2))

}
}
