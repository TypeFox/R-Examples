#' Probability density function (pdf) for an uncensored epanechnikov distribution
#'
#' This function is simply a polynomial of second degree.
#'
#' @param x point on x-axis.
#' @param mu mean of distribution.
#' @param r half the range of the distribution, ie the distance from the mean to the smallest/largest value supported by the distribution. \code{r=5^.5} corresponds to a standard deviation of 1.
#' @return point density associated with \code{x}, \code{mu} and \code{r}.
#' @keywords distribution
#' @examples
#' #Probability distribution function, epanechnikov:
#' curve(depan(x),col='blue',ylim=c(0,.4),xlim=c(-3.5,3.5),yaxs='i',xaxs='i',
#' main='Probability distribution function',ylab='Probability')
#'
#' #Probability distribution function, normal:
#' curve(dnorm(x),col='green',add=TRUE)
#'
#' #Legend
#' legend(x=-3.5,y=.4,legend=c('Epanechnikov pdf','Normal pdf'),lty=c(1,1),col=c('blue','green'))


depan <- function(x = 0, mu = 0, r = 5^0.5) {
    # Distribution function
    if (any(r <= 0)) {
        stop("Range must be strictly positive")
    }
    ifelse(abs(x - mu) < r, 3/4 * (1 - ((x - mu)/r)^2)/r, 0)
    
} 
