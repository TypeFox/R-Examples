#' Generate random uncensored epanechnikov-distributed data
#'
#' This function works in conjuncture with \code{qepan} and \code{runif}
#'
#' @param n number of data points.
#' @param mu mean of distribution.
#' @param r half the range of the distribution, ie the distance from the mean to the smallest/largest value supported by the distribution. \code{r=5^.5} corresponds to a standard deviation of 1.
#' @return vector of random variables.
#' @keywords distribution
#' @examples
#' #Generate and plot 10000 random observations:
#' hist(repan(10000,mu=100,r=10))


repan <- function(n, mu = 0, r = 5^0.5) {
    if (any(r <= 0)) {
        stop("Range must be strictly positive")
    }
    if (n%%1 != 0 | n < 1) {
        stop("n must be an integer")
    }
    
    qepan(runif(n, min = 0, max = 1), mu, r)
} 
