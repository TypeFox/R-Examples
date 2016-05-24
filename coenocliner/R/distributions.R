##' @title Wrappers to random number generators for use with coenocliner
##'
##' @description These functions are simple wrappers around existing random number generators in R to provide stochasic count data for simulated species.
##'
##' @references Bolker, B.M. (2008) \emph{Ecological Models and Data
##' in R.} Princeton University Press.
##'
##' @param n the number of random draws, equal to number of species times the number of gradient locations.
##' @param mu the mean or expectation of the distribution. For \code{Bernoulli}, \code{Binomial}, and \code{BetaBinomial()} this is the probability of occurrence as given by the response function.
##' @param alpha numeric; parameter for the negative binomial distribution.
##'
##' @return a vector of random draws from the stated distribution.
##'
##' @author Gavin L. Simpson
##'
##' @rdname distributions
##'
##' @name distributions
##'
##' @keywords distribution
##'
##' @importFrom stats rpois rgamma
`NegBin` <- function(n, mu, alpha) {
    mu <- mu * rgamma(n, shape = alpha, rate = 1/alpha)
    rpois(n, lambda = mu)
}

##' @rdname distributions
##'
##' @importFrom stats rpois
`Poisson` <- function(n, mu) {
    rpois(n, lambda = mu)
}

##' @rdname distributions
##'
##' @importFrom stats rbinom
`Bernoulli` <- function(n, mu) {
    rbinom(n, size = 1, prob = mu)
}

##' @rdname distributions
##'
##' @importFrom stats rbinom
##'
##' @param size numeric; binomial denominator, the total number of individuals counted for example
`Binomial` <- function(n, mu, size) {
    rbinom(n = n, size = size, prob = mu)
}

##' @rdname distributions
##'
##' @param theta numeric; a positive \emph{inverse} overdispersion parameter for the Beta-Binomial distribution. Low values give high overdispersion. The variance is  \code{size*mu*(1-mu)*(1+(size-1)/(theta+1))} (Bolker, 2008)
##'
##' @importFrom stats rbeta rbinom
`BetaBinomial` <- function(n, mu, size, theta) {
    ## follows Bolker (2008) and derive a and b from mu and theta
    ## mu == pi (or prob, or p)
    a <- theta * mu
    b <- theta * (1 - mu)
    rbinom(n = n, size = size,
           prob = rbeta(n = n, shape1 = a, shape2 = b))
}

##' @rdname distributions
##'
##' @param zprobs numeric; zero-inflation parameter giving the proportion of extraneous zeros. Must be in range \eqn{0 \dots 1}{0 to 1}.
##' @importFrom stats runif rpois
`ZIP` <- function(n, mu, zprobs) {
    ifelse(runif(n) > zprobs, rpois(n, lambda = mu), 0)
}

##' @rdname distributions
##'
##' @importFrom stats rpois rgamma runif
`ZINB` <- function(n, mu, alpha, zprobs) {
    ifelse(runif(n) > zprobs,
           rpois(n, lambda = mu *
                 rgamma(n, shape = alpha, rate = 1/alpha)),
           0)
}

##' @rdname distributions
##'
##' @importFrom stats rbinom runif
## Zero-inflated Binomial
`ZIB` <- function(n, mu, size, zprobs) {
    ifelse(runif(n) > zprobs,
           rbinom(n, size = size, prob = mu),
           0)
}

##' @rdname distributions
##'
##' @importFrom stats runif
## Zero-inflated Beta-Binomial
`ZIBB` <- function(n, mu, size, theta, zprobs) {
    ifelse(runif(n) > zprobs,
           BetaBinomial(n, mu = mu, size = size, theta = theta),
           0)
}
