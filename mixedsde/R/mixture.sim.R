#' Simulation Of A Mixture Of Two Normal Or Gamma Distributions
#' 
#' @description Simulation of M random variables from a mixture of two Gaussian or Gamma distributions
#' @param M number of simulated variables 
#' @param density.phi name of the chosen density 'mixture.normal' or 'mixture.gamma'
#' @param param vector of parameters with the proportion of mixture of the two distributions and means and standard-deviations of the two normal or 
#'   shapes and scales of the two Gamma distribution 
#' @return
#' \item{Y}{vector of simulated variables}
#' @details
#' If 'mixture.normal', the distribution is \eqn{p N(\mu1,\sigma1^2) + (1-p)N(\mu2, \sigma2^2)}
#' 
#' and param=c(p, \eqn{\mu1, \sigma1, \mu2, \sigma2})
#' 
#' If 'mixture.gamma', the distribution is \eqn{p Gamma(shape1,scale1) + (1-p)Gamma(shape2,scale2)}
#' 
#' and param=c(p, shape1, scale1, shape2, scale2)
#' @examples 
#' density.phi <- 'mixture.gamma'
#' param <- c(0.2,1.8,0.5,5.05,1); M <- 200
#' gridf <- seq(0, 8, length = 200)  
#' f <- param[1] * 1/gamma(param[2]) * (gridf)^(param[2]-1) * 
#'            exp(-(gridf) / param[3]) / param[3]^param[2] + 
#'	(1-param[1]) * 1/gamma(param[4]) * (gridf)^(param[4]-1) * 
#'	    exp(-(gridf) / param[5]) / param[5]^param[4]
#' Y <- mixture.sim(M, density.phi, param)
#' hist(Y)
#' lines(gridf, f)


mixture.sim <- function(M, density.phi, param) {
    
    propmixed <- param[1]
    p1 <- rbinom(M, 1, propmixed)
    
    if (density.phi == "mixture.normal") {
        mu1 <- param[2]
        omega1 <- param[3]
        mu2 <- param[4]
        omega2 <- param[5]
        Y <- p1 * rnorm(M, mean = mu1, sd = omega1) + (1 - p1) * rnorm(M, mean = mu2, sd = omega2)
    }
    if (density.phi == "mixture.gamma") {
        shape1 <- param[2]
        scale1 <- param[3]
        shape2 <- param[4]
        scale2 <- param[5]
        Y <- p1 * rgamma(M, shape = shape1, scale = scale1) + (1 - p1) * rgamma(M, shape = shape2, scale = scale2)
    }
    
    return(Y)
} 
