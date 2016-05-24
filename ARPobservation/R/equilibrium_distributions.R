#' @title Constructor for class \code{eq_dist}
#' 
#' @description
#' The \code{eq_dist} class consists of a pair of component functions for generating random variates from a 
#' specified distribution and the corresponding equilibrium distribution.
#' 
#' @param r_gen function for generating random deviates.
#' @param r_eq function for generating random deviates from the corresponding equilibrium distribution.
#' 
#' @details Both functions must take arguments \code{n} and \code{mean}. Currently, the following distributions
#' are implemented:
#' \itemize{
#' \item \code{\link{F_exp}} - Exponential
#' \item \code{\link{F_gam}} - Gamma
#' \item \code{\link{F_gam_mix}} - Mixture of two gammas
#' \item \code{\link{F_weib}} - Weibull
#' \item \code{\link{F_unif}} - Uniform
#' \item \code{\link{F_const}} - Constant
#' }
#' 
#' @return Object of class \code{eq_dist} with components \code{r_gen} and \code{r_eq}.
#' 
#' @export 

eq_dist <- function(r_gen, r_eq) {
  eq_dist <- list(r_gen = r_gen, r_eq = r_eq)
  class(eq_dist) <- "eq_dist"
  return(eq_dist)
}


## exponential distribution ####

#' @title Exponential distribution and related equilibrium distribution
#' 
#' @description
#' Random number generation from exponential distributions, for use with \code{\link{r_behavior_stream}}.
#' 
#' @return Object of class \code{\link{eq_dist}} with components \code{r_gen} and \code{r_eq}.
#' 
#' The function \code{r_gen(n, mean)} generates random deviates from an exponential distribution with specified \code{mean}.
#' 
#' The function \code{r_eq(n, mean)} generates random deviates from an exponential distribution with specified \code{mean}.

#' 
#' @examples
#' hist(F_exp()$r_gen(1000, 3))
#' hist(F_exp()$r_eq(1000, 3))
#' 
#' @export 

F_exp <- function()
  eq_dist(r_gen = function(n, mean) rexp(n, rate = 1 / mean),
          r_eq = function(n, mean) rexp(n, rate = 1 / mean))





## gamma distribution ####

pgamma_eq <- function(x, mean, shape) 
  x / mean + pgamma(x, shape = shape + 1, scale = mean / shape) - 
  pgamma(x, shape = shape, scale = mean / shape) * x / mean

rgamma_eq <- function(n, mean, shape) mapply(function(p, m) 
  uniroot(function(y) p - pgamma_eq(y, mean = m, shape = shape), lower = 0, upper = m * 10^5)$root,
                                             p = runif(n), m = mean)

#' @title Gamma distribution and related equilibrium distribution
#' 
#' @description
#' Random number generation from a gamma distribution and the related equilibrium distribution, 
#' for use with \code{\link{r_behavior_stream}}.
#' 
#' @param shape shape parameter
#' 
#' @return Object of class \code{\link{eq_dist}} with components \code{r_gen} and \code{r_eq}.
#' 
#' The function \code{r_gen(n, mean)} generates random deviates from a gamma distribution with specified 
#' \code{mean} and \code{shape} parameters.
#' 
#' The function \code{r_eq(n, mean)} generates random deviates from the equilibrium distribuion corresponding
#' to the gamma distribution with specified \code{mean} and \code{shape} parameters.
#' 
#' @examples
#' hist(F_gam(2)$r_gen(1000, 3))
#' hist(F_gam(2)$r_eq(1000, 3))
#' 
#' @export 

F_gam <- function(shape)
  eq_dist(r_gen = function(n, mean) rgamma(n, shape = shape, scale = mean / shape),
          r_eq = function(n, mean) rgamma_eq(n, mean, shape))





## gamma mixture ####

pgamma_mix_eq <- function(x, mean, shape1, shape2, scale_ratio, mix) {
  theta2 <- mean / (mix * shape1 * scale_ratio + (1 - mix) * shape2)
  theta1 <- scale_ratio * theta2
  (mix * shape1 * theta1 * pgamma_eq(x, shape1 * theta1, shape1) + 
     (1 - mix) * shape2 * theta2 * pgamma_eq(x, shape2 * theta2, shape2)) / mean
}

rgamma_mix_eq <- function(n, mean, shape1, shape2, scale_ratio, mix) 
  mapply(function(p, m)
    uniroot(function(y) p - pgamma_mix_eq(y, mean=m, shape1, shape2, scale_ratio, mix), 
          lower = 0, upper = m * 10^5)$root,
          p = runif(n), m = mean)

#' @title Mixture of two gamma distributions and related equilibrium distribution
#' 
#' @description
#' Random number generation from a mixture of two gamma distributions and the related equilibrium distribution, 
#' for use with \code{\link{r_behavior_stream}}.
#' 
#' @param shape1 shape parameter for first mixture component, \eqn{k_1}
#' @param shape2 shape parameter for second mixture component, \eqn{k_2}
#' @param scale_ratio ratio of first scale component to second scale component, \eqn{\theta_1 / \theta_2}
#' @param mix mixing proportion of first component, \eqn{p}
#' 
#' @return Object of class \code{\link{eq_dist}} with components \code{r_gen} and \code{r_eq}.
#' 
#' The function \code{r_gen(n, mean)} generates random deviates from a mixture of two gamma distributions with specified 
#' \code{mean}, \code{shape1}, \code{shape2}, \code{scale_ratio}, and \code{mix}. The cumulative distribution function 
#' is given by \deqn{F(x) = p \Gamma(x; k_1, \theta_1) + (1 - p) \Gamma(x; k_2, \theta_2),} where \eqn{\Gamma(x; k, \theta)}
#' is the cumulative distribution function of a Gamma random variable with shape \eqn{k} and scale \eqn{\theta}, and
#' the scale parameters are determined by the specified \code{mean} and \code{scale_ratio}.
#' 
#' The function \code{r_eq(n, mean)} generates random deviates from the equilibrium distribuion corresponding
#' to the mixture of gamma distributions.
#' 
#' @examples
#' hist(F_gam_mix(2, 2, 1 / 12, 3 / 5)$r_gen(1000, 20))
#' hist(F_gam_mix(2, 2, 1 / 12, 3 / 5)$r_eq(1000, 20))
#' 
#' @export 

F_gam_mix <- function(shape1, shape2, scale_ratio, mix)
  eq_dist(r_gen = function(n, mean) {
    m <- rbinom(n, 1, mix)
    shape <- c(shape1, shape2)[2 - m]
    scale <- c(scale_ratio, 1)[2 - m] * mean / (mix * shape1 * scale_ratio + (1 - mix) * shape2) 
    rgamma(n, shape=shape, scale=scale)
  }, 
  r_eq = function(n, mean) rgamma_mix_eq(n, mean, shape1, shape2, scale_ratio, mix))


## Weibull distribution ####

pweibull_eq <- function(x, mean, shape) {
  scale <- mean / gamma(1 + 1 / shape)
  integrate(function(z) exp(-(z / scale)^shape), 0, x)$value / mean
}

rweibull_eq <- function(n, mean, shape) mapply(function(p, m) 
  uniroot(function(y) p - pweibull_eq(y, mean = m, shape = shape), lower = 0, upper = m * 10^3)$root,
  p = runif(n), m = mean)

#' @title Weibull distribution and related equilibrium distribution
#' 
#' @description
#' Random number generation from a Weibull distribution and the related equilibrium distribution, 
#' for use with \code{\link{r_behavior_stream}}.
#' 
#' @param shape shape parameter
#' 
#' @return Object of class \code{\link{eq_dist}} with components \code{r_gen} and \code{r_eq}.
#' 
#' The function \code{r_gen(n, mean)} generates random deviates from a Weibull distribution with specified 
#' \code{mean} and \code{shape} parameters.
#' 
#' The function \code{r_eq(n, mean)} generates random deviates from the equilibrium distribuion corresponding
#' to the Weibull distribution with specified \code{mean} and \code{shape} parameters.
#' 
#' @examples
#' hist(F_gam(2)$r_gen(1000, 3))
#' hist(F_gam(2)$r_eq(1000, 3))
#' 
#' @export 

F_weib <- function(shape) 
  eq_dist(r_gen = function(n, mean) rweibull(n, shape = shape, scale = mean / gamma(1 + 1 / shape)),
          r_eq = function(n, mean) rweibull_eq(n, mean, shape))



## uniform distribution on (0, 2 * mean) ####

#' @title Uniform distribution and related equilibrium distribution
#' 
#' @description
#' Random number generation from a uniform distribution and the related equilibrium distribution, 
#' for use with \code{\link{r_behavior_stream}}.
#' 
#' @return Object of class \code{\link{eq_dist}} with components \code{r_gen} and \code{r_eq}.
#' 
#' The function \code{r_gen(n, mean)} generates random deviates from a uniform distribution with specified 
#' \code{mean} \eqn{\mu} on the interval \eqn{(0, 2 \mu)}. The cumulative distribution function 
#' is given by \eqn{F(x) = x / 2 \mu}.
#' 
#' The function \code{r_eq(n, mean)} generates random deviates from the equilibrium distribuion corresponding
#' to a uniform distribution on the interval \eqn{(0, 2 \mu)}. The cumulative distribution function is given by 
#' \deqn{F(x) = x (4 \mu - x) / (4 \mu^2).}
#' 
#' @examples
#' hist(F_unif()$r_gen(1000, 2))
#' hist(F_unif()$r_eq(1000, 2))
#' 
#' @export 

F_unif <- function()
  eq_dist(r_gen = function(n, mean) runif(n, min = 0, max = 2 * mean),
          r_eq = function(n, mean) 2 * mean * (1 - sqrt(1 - runif(n))))





## constant (degenerate) distribution ####

#' @title Constant (degenerate) distribution and related equilibrium distribution
#' 
#' @description
#' Generation from a degenerate distribution and random number generation from the related equilibrium distribution, 
#' for use with \code{\link{r_behavior_stream}}.
#' 
#' @return Object of class \code{\link{eq_dist}} with components \code{r_gen} and \code{r_eq}.
#' 
#' The function \code{r_gen(n, mean)} simply returns a vector of length \code{n} with all values equal to \code{mean}.
#' 
#' The function \code{r_eq(n, mean)} generates random deviates from a uniform distribution on the interval (0, mean).
#' 
#' @examples
#' hist(F_const()$r_gen(1000, 2))
#' hist(F_const()$r_eq(1000, 2))
#' 
#' @export 

F_const <- function() 
  eq_dist(r_gen = function(n, mean) rep(mean, length.out = n),
          r_eq = function(n, mean) runif(n, min=0, max=mean))
  

