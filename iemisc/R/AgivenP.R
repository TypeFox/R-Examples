#' Annual value given Present value (Engineering Economics)
#'
#' Compute A given P
#'
#' A is expressed as
#'
#' 	\deqn{A = P\left[\frac{i\left(1 + i\right)^n}{\left(1 + i\right)^n - 1}\right]}
#'
#' \describe{
#'	\item{\emph{A}}{the "uniform series amount (occurs at the end of each
#'     interest period)"}
#'	\item{\emph{P}}{the "present equivalent"}
#'	\item{\emph{i}}{the "effective interest rate per interest period"}
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
#' @return AgivenP numeric vector that contains the annual value(s) rounded
#'    to 2 decimal places
#' @return AP data.frame of both n (0 to n) and the resulting annual values
#'    rounded to 2 decimal places
#'
#' @references
#' William G. Sullivan, Elin M. Wicks, and C. Patrick Koelling, \emph{Engineering Economy}, Fourteenth Edition, Upper Saddle River, New Jersey: Pearson/Prentice Hall, 2009, page 136, 142, 164, 166.
#'
#' @encoding UTF-8
#'
#'
#' @examples
#' library(iemisc)
#' # Example for equation 4-14 from the Reference text (page 136)
#' AgivenP(17000, 4, 1, "annual")
#' # the interest rate is 1% per month and n is 4 months
#'
#' AP(17000, 4, 1, "annual")
#' # the interest rate is 1% per month and n is 4 months
#'
#'
#' # Example 4-30 from the Reference text (page 166)
#' AgivenP(10000, 5, 12, "month")
#' # the interest rate is 12% compounded monthly for 5 years
#'
#' AP(10000, 5, 12, "month")
#' # the interest rate is 12% compounded monthly for 5 years
#'
#'
#'
#' @import data.table
#'
#' @name AgivenP
NULL

#' @export
#' @rdname AgivenP
AgivenP <- function (P, n, i, frequency = c("annual", "semiannual", "quarter", "bimonth", "month", "daily")) {

i <- i / 100

fr <- frequency

if (fr == "annual") {
fr <- 1
n <- n * fr

i <- i / fr

AgivenP <- P * ((i * ((1 + i) ^ n)) / (((1 + i) ^ n) - 1))

return(round(AgivenP, digits = 2))

} else if (fr == "semiannual") {

fr <- 2
n <- n * fr

i <- i / fr

AgivenP <- P * ((i * ((1 + i) ^ n)) / (((1 + i) ^ n) - 1))

return(round(AgivenP, digits = 2))

} else if (fr == "quarter") {

fr <- 4
n <- n * fr

i <- i / fr

AgivenP <- P * ((i * ((1 + i) ^ n)) / (((1 + i) ^ n) - 1))

return(round(AgivenP, digits = 2))

} else if (fr == "bimonth") {

fr <- 6
n <- n * fr

i <- i / fr

AgivenP <- P * ((i * ((1 + i) ^ n)) / (((1 + i) ^ n) - 1))

return(round(AgivenP, digits = 2))

} else if (fr == "month") {

fr <- 12
n <- n * fr

i <- i / fr

AgivenP <- P * ((i * ((1 + i) ^ n)) / (((1 + i) ^ n) - 1))

return(round(AgivenP, digits = 2))

} else if (fr == "daily") {

fr <- 365
n <- n * fr

i <- i / fr

AgivenP <- P * ((i * ((1 + i) ^ n)) / (((1 + i) ^ n) - 1))

return(round(AgivenP, digits = 2))

}
}



#' @export
#' @rdname AgivenP
AP <- function (P, n, i, frequency = c("annual", "semiannual", "quarter", "bimonth", "month", "daily")) {

i <- i / 100

fr <- frequency

if (fr == "annual") {
fr <- 1
n <- n * fr

AP <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
AP[[y]] <- P * ((i * ((1 + i) ^ seq(n))) / (((1 + i) ^ seq(n)) - 1))
}

AP <- data.table(seq(n), unlist(AP))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

AP <- rbind(P0, AP)

setnames(AP, c("n (periods)", "Annual Worth ($US)"))
AP <- setDF(AP)

return(round(AP, digits = 2))

} else if (fr == "semiannual") {

fr <- 2
n <- n * fr

i <- i / fr

AP <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
AP[[y]] <- P * ((i * ((1 + i) ^ seq(n))) / (((1 + i) ^ seq(n)) - 1))
}

AP <- data.table(seq(n), unlist(AP))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

AP <- rbind(P0, AP)

setnames(AP, c("n (periods)", "Annual Worth ($US)"))
AP <- setDF(AP)

return(round(AP, digits = 2))

} else if (fr == "quarter") {

fr <- 4
n <- n * fr

i <- i / fr

AP <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
AP[[y]] <- P * ((i * ((1 + i) ^ seq(n))) / (((1 + i) ^ seq(n)) - 1))
}

AP <- data.table(seq(n), unlist(AP))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

AP <- rbind(P0, AP)

setnames(AP, c("n (periods)", "Annual Worth ($US)"))
AP <- setDF(AP)

return(round(AP, digits = 2))

} else if (fr == "bimonth") {

fr <- 6
n <- n * fr

i <- i / fr

AP <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
AP[[y]] <- P * ((i * ((1 + i) ^ seq(n))) / (((1 + i) ^ seq(n)) - 1))
}

AP <- data.table(seq(n), unlist(AP))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

AP <- rbind(P0, AP)

setnames(AP, c("n (periods)", "Annual Worth ($US)"))
AP <- setDF(AP)

return(round(AP, digits = 2))

} else if (fr == "month") {

fr <- 12
n <- n * fr

i <- i / fr

AP <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
AP[[y]] <- P * ((i * ((1 + i) ^ seq(n))) / (((1 + i) ^ seq(n)) - 1))
}

AP <- data.table(seq(n), unlist(AP))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

AP <- rbind(P0, AP)

setnames(AP, c("n (periods)", "Annual Worth ($US)"))
AP <- setDF(AP)

return(round(AP, digits = 2))

} else if (fr == "daily") {

fr <- 365
n <- n * fr

i <- i / fr

AP <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
AP[[y]] <- P * ((i * ((1 + i) ^ seq(n))) / (((1 + i) ^ seq(n)) - 1))
}

AP <- data.table(seq(n), unlist(AP))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

AP <- rbind(P0, AP)

setnames(AP, c("n (periods)", "Annual Worth ($US)"))
AP <- setDF(AP)

return(round(AP, digits = 2))

}
}
