#' pSnbinom returns the distribution function for the sum of random variable with negative binomial distributions
#' @title Distribution function for the sum of random variable with negative binomial distributions. 
#' @author Marc Girondot
#' @return pSnbinom returns distribution function
#' @param q vector of quantiles.
#' @param size target for number of successful trials, or dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.
#' @param prob probability of success in each trial. 0 < prob <= 1.
#' @param mu alternative parametrization via mean.
#' @param log.p	logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @param infinite Number of maximal iterations; check different values to determine the error in estimation.
#' @description Distribution function for the sum of random variable with negative binomial distributions.
#' @family Distribution of sum of random variable with negative binomial distributions
#' @examples
#' alpha <- c(2.1, 2.05, 2)
#' mu <- c(10, 30, 20)
#' p <- pSnbinom(q=10, size=alpha, mu=mu, lower.tail = TRUE)
#' @export

pSnbinom <- function(q=stop("At least one quantile must be provided"), 
                     size=NULL, 
                     prob=NULL, mu=NULL, lower.tail = TRUE, log.p = FALSE, infinite=100) {
  
  # prob=NULL; mu=NULL; log = FALSE; infinite=10
  
  if (is.null(mu) + is.null(size) + is.null(prob) != 1) stop("Two values among mu, size and prob must be provided")
  
  m <- max(c(length(size), length(prob), length(mu)))
  if (!is.null(mu)) mu <- rep(mu, m)[1:m]
  if (!is.null(size)) size <- rep(size, m)[1:m]
  if (!is.null(prob)) prob <- rep(prob, m)[1:m]
  
  if (is.null(prob)) prob <- size/(size+mu)
  if (is.null(mu)) mu <- size/prob - size
  if (is.null(size))  size  <- (prob * mu) / (1 - prob)
  
#  if (length(prob)<length(size)) prob <- rep(prob, length(size))[1:length(size)]
#  if (length(size)<length(prob)) size <- rep(size, length(prob))[1:length(prob)]
  
  pp <- vapply(q, FUN=function(qq) {
  
  l <- dSnbinom(0:qq, prob=prob, size=size, mu=NULL, log=FALSE, infinite=infinite)
  p <- sum(l)
  if (!lower.tail) p <- 1-p
  if (log.p) p <- log(p)
  
  return(p)
  }, FUN.VALUE = 10.1)
  
  return(pp)
  }
  