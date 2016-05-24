#' Calculate expected value of censored variable
#'
#' This function calculates the expected value of a random censored epanechnikov-distributed variable with a given censoring point.
#' The inverse of this function is \code{cepan}.
#'
#' @param c censoring point.
#' @param mu mean of distribution prior to censoring.
#' @param r half the range of the distribution, ie the distance from the mean to the smallest/largest value supported by the distribution. \code{r=5^.5} corresponds to a standard deviation of 1.
#' @param side_censored indicates whether the variable is \code{left} or \code{right} censored. Default is \code{side_censored='left'}
#' @return the expected value associated with \code{c}, \code{mu} and \code{r}.
#' @keywords distribution
#' @examples
#' #Expected value of an epan-distributed variable left-censored at 100 (given mu=100 and r=10):
#' evepan(c=100,mu=100,r=10)
#'
#' #Expected value as a function of censoring point, epanechnikov distribution:
#' curve(evepan(c=x),col='blue',xlim=c(-sqrt(5),sqrt(5)),yaxs='i',xaxs='i',
#' main='Expected value as a function of censoring point',xlab='Censoring point',ylab='Expected value')
#'
#' #Expected value as a function of censoring point, normal distribution:
#' curve(dnorm(x)+pnorm(x)*x,col='green',add=TRUE)
#'
#' #Expected value as a function of censoring point, no uncertainty:
#' curve(1*x,col='grey',add=TRUE)
#'
#' #Legend
#' legend(x=-sqrt(5),y=sqrt(5),legend=c('Epanechnikov','Normal distribution','No uncertainty'),
#' lty=c(1,1),col=c('blue','green','grey'))


evepan <- function(c = 0, mu = 0, r = 5^0.5, side_censored = "left") {
    # Expected abatement
    if (any(r <= 0)) {
        stop("Range must be strictly positive")
    }

    if (any(!(side_censored %in% c("left", "right")))) {
        stop("side_censored must either be 'left' or 'right'")
    }

    setsign <- ifelse(side_censored == "left", 1, -1)

    alpha <- setsign * (c - mu)/r

    ev <- setsign * r/16 * (1 - alpha)^3 * (3 + alpha) + c

    ifelse(abs(alpha) <= 1, ev, ifelse(alpha < -1, mu, c))

}
