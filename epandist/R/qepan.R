#' Quantile function for an uncensored epanechnikov distribution
#'
#' The inverse of this function is \code{pepan}.
#'
#' @param p probability.
#' @param mu mean of distribution.
#' @param r half the range of the distribution, ie the distance from the mean to the smallest/largest value supported by the distribution. \code{r=5^.5} corresponds to a standard deviation of 1.
#' @return the quantile associated with \code{x}, \code{mu} and \code{r}.
#' @keywords distribution
#' @examples
#' #Calculate the lower quartile of an epan-distributed variable:
#' qepan(p=.25,mu=0,r=sqrt(5))
#' 
#' #Use qepan to confirm analytical solution
#' #Find the quantile corresponding to p=(5+sqrt(5))/8=.9045 when mu=0 and r=sqrt(5):
#' qepan(p=(5+sqrt(5))/8,mu=0,r=sqrt(5))
#' #This is equal to 
#' (5-sqrt(5))/2 


qepan <- function(p, mu = 0, r = 5^0.5) {
    if (any(r <= 0)) {
        stop("Range must be strictly positive")
    }
    if (any(abs(p - 0.5) > 0.5)) {
        stop("p must be between 0 and 1")
    }
    
    2 * sin(asin(2 * p - 1)/3) * r + mu
} 
