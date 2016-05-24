#' Present value given Annual value (Engineering Economics)
#'
#' Compute P given A
#'
#' P is expressed as
#'
#' 	\deqn{P = A\left[\frac{\left(1 + i\right)^n - 1}{i\left(1 + i\right)^n}\right]}
#'
#' \describe{
#'	\item{\emph{P}}{the "present equivalent"}
#'	\item{\emph{A}}{the "uniform series amount (occurs at the end of each
#'     interest period)"}
#'	\item{\emph{i}}{the "effective interest rate per interest period"}
#'	\item{\emph{n}}{the "number of interest periods"}
#' }
#'
#'
#' @param A numeric vector that contains the annual value(s)
#' @param n numeric vector that contains the period value(s)
#' @param i numeric vector that contains the interest rate(s) as a percent
#' @param frequency character vector that contains the frequency used to
#'    obtain the number of periods [annual (1), semiannual (2), quarter (4),
#'    bimonth (6), month (12), daily (365)]
#'
#' @return PgivenA numeric vector that contains the present value(s) rounded
#'    to 2 decimal places
#' @return PA data.frame of both n (0 to n) and the resulting present values
#'    rounded to 2 decimal places
#'
#' @references
#' William G. Sullivan, Elin M. Wicks, and C. Patrick Koelling, \emph{Engineering Economy}, Fourteenth Edition, Upper Saddle River, New Jersey: Pearson/Prentice Hall, 2009, page 133-134, 142, 164.
#'
#' @encoding UTF-8
#'
#'
#' @examples
#' library(iemisc)
#' # Example 4-9 from the Reference text (page 133-134)
#' PgivenA(20000, 5, 15, "annual") # the interest rate is 15%
#'
#' PA(20000, 5, 15, "annual") # the interest rate is 15%
#'
#'
#' @import data.table
#'
#' @name PgivenA
NULL

#' @export
#' @rdname PgivenA
PgivenA <- function (A, n, i, frequency = c("annual", "semiannual", "quarter", "bimonth", "month", "daily")) {

i <- i / 100

fr <- frequency

if (fr == "annual") {
fr <- 1
n <- n * fr

i <- i / fr

PgivenA <- A * (((1 + i) ^ n - 1) / (i * ((1 + i) ^ n)))

return(round(PgivenA, digits = 2))

} else if (fr == "semiannual") {

fr <- 2
n <- n * fr

i <- i / fr

PgivenA <- A * (((1 + i) ^ n - 1) / (i * ((1 + i) ^ n)))

return(round(PgivenA, digits = 2))

} else if (fr == "quarter") {

fr <- 4
n <- n * fr

i <- i / fr

PgivenA <- A * (((1 + i) ^ n - 1) / (i * ((1 + i) ^ n)))

return(round(PgivenA, digits = 2))

} else if (fr == "bimonth") {

fr <- 6
n <- n * fr

i <- i / fr

PgivenA <- A * (((1 + i) ^ n - 1) / (i * ((1 + i) ^ n)))

return(round(PgivenA, digits = 2))

} else if (fr == "month") {

fr <- 12
n <- n * fr

i <- i / fr

PgivenA <- A * (((1 + i) ^ n - 1) / (i * ((1 + i) ^ n)))

return(round(PgivenA, digits = 2))

} else if (fr == "daily") {

fr <- 365
n <- n * fr

i <- i / fr

PgivenA <- A * (((1 + i) ^ n - 1) / (i * ((1 + i) ^ n)))

return(round(PgivenA, digits = 2))

}
}



#' @export
#' @rdname PgivenA
PA <- function (A, n, i, frequency = c("annual", "semiannual", "quarter", "bimonth", "month", "daily")) {

i <- i / 100

fr <- frequency

if (fr == "annual") {
fr <- 1
n <- n * fr

PA <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
PA[[y]] <- A * (((1 + i) ^ seq(n) - 1) / (i * ((1 + i) ^ seq(n))))
}

PA <- data.table(seq(n), unlist(PA))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

PA <- rbind(P0, PA)

setnames(PA, c("n (periods)", "Present Worth ($US)"))
PA <- setDF(PA)

return(round(PA, digits = 2))

} else if (fr == "semiannual") {

fr <- 2
n <- n * fr

i <- i / fr

PA <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
PA[[y]] <- A * (((1 + i) ^ seq(n) - 1) / (i * ((1 + i) ^ seq(n))))
}

PA <- data.table(seq(n), unlist(PA))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

PA <- rbind(P0, PA)

setnames(PA, c("n (periods)", "Present Worth ($US)"))
PA <- setDF(PA)

return(round(PA, digits = 2))

} else if (fr == "quarter") {

fr <- 4
n <- n * fr

i <- i / fr

PA <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
PA[[y]] <- A * (((1 + i) ^ seq(n) - 1) / (i * ((1 + i) ^ seq(n))))
}

PA <- data.table(seq(n), unlist(PA))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

PA <- rbind(P0, PA)

setnames(PA, c("n (periods)", "Present Worth ($US)"))
PA <- setDF(PA)

return(round(PA, digits = 2))

} else if (fr == "bimonth") {

fr <- 6
n <- n * fr

i <- i / fr

PA <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
PA[[y]] <- A * (((1 + i) ^ seq(n) - 1) / (i * ((1 + i) ^ seq(n))))
}

PA <- data.table(seq(n), unlist(PA))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

PA <- rbind(P0, PA)

setnames(PA, c("n (periods)", "Present Worth ($US)"))
PA <- setDF(PA)

return(round(PA, digits = 2))

} else if (fr == "month") {

fr <- 12
n <- n * fr

i <- i / fr

PA <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
PA[[y]] <- A * (((1 + i) ^ seq(n) - 1) / (i * ((1 + i) ^ seq(n))))
}

PA <- data.table(seq(n), unlist(PA))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

PA <- rbind(P0, PA)

setnames(PA, c("n (periods)", "Present Worth ($US)"))
PA <- setDF(PA)

return(round(PA, digits = 2))

} else if (fr == "daily") {

fr <- 365
n <- n * fr

i <- i / fr

PA <- vector("list", length(1:n))
## Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (y in 1:length(n)) {
PA[[y]] <- A * (((1 + i) ^ seq(n) - 1) / (i * ((1 + i) ^ seq(n))))
}

PA <- data.table(seq(n), unlist(PA))

P0 <- NA
P0 <- data.table(0, P0)
P0 <- setnames(P0, 2, "V2")

PA <- rbind(P0, PA)

setnames(PA, c("n (periods)", "Present Worth ($US)"))
PA <- setDF(PA)

return(round(PA, digits = 2))

}
}
