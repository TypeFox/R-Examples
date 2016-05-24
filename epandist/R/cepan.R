#' Calculate censoring point
#'
#'This function calculates the censoring point of a random censored epanechnikov-distributed variable associated a given expected value.
#'The inverse of this function is \code{evepan}.
#'
#' @param ev expected value.
#' @param mu mean of distribution prior to censoring.
#' @param r half the range of the distribution, ie the distance from the mean to the smallest/largest value supported by the distribution. \code{r=5^.5} corresponds to a standard deviation of 1.
#' @param side_censored indicates whether the variable is \code{left} or \code{right} censored. Default is \code{side_censored='left'}
#' @return the censoring point associated with \code{ev}, \code{mu} and \code{r}.
#' @keywords distribution
#' @examples
#' #Censoring point of a left-censored epan-distributed variable
#' #with an expected value of 3 (given mu=0 and r=16):
#' cepan(ev=3,mu=0,r=16)
#'
#' #Censoring point of a right-censored epan-distributed variable
#' #with an expected value of 103 (given mu=100 and r=32):
#' cepan(ev=94,mu=100,r=32,side_censored="right")
#' #Results are usually not an integer though and rarely coinciding with mu




cepan <- function(ev, mu = 0, r = 5^0.5, side_censored = "left") {

    if (any(r <= 0)) {
        stop("Range must be strictly positive")
    }

    if (any(!(side_censored %in% c("left", "right")))) {
        stop("side_censored must either be 'left' or 'right'")
    }

    setsign <- ifelse(side_censored == "left", 1, -1)


    if (any(setsign * ev <= setsign * mu)) {
        stop("The expected value of a left-censored variable must be greater than the mean of the uncensored distribution.\n         Likewise the expected value of a right-censored variable must be less than the mean of the uncensored distribution.")
    }


    k <- setsign * (ev - mu)/r

    tempvar0 <- (3^(1/2) * (27 * k^2 - 16 * k^3)^(1/2) + 9 * k)^(1/3)
    tempvar1 <- (2 * 2^(2/3) * k/3^(1/3)/tempvar0 + 1 * 2^(1/3) * tempvar0/3^(2/3) + 1)^(1/2)
    tempvar2 <- (-8 * 2^(2/3) * k/3^(1/3)/tempvar0 - 4 * 2^(1/3) * tempvar0/3^(2/3) + 8/tempvar1 + 8)^(1/2)/2
    solution_to_quad <- tempvar1 - tempvar2

    solution_in_interval <- setsign * solution_to_quad * r + mu

    ifelse(setsign * (ev - mu) < r, solution_in_interval, ev)


}
