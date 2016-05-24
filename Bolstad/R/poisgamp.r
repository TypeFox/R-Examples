#' Poisson sampling with a gamma prior
#' 
#' Evaluates and plots the posterior density for \eqn{\mu}{mu}, the mean rate
#' of occurance in a Poisson process and a \eqn{gamma} prior on \eqn{\mu}{mu}
#' 
#' 
#' @param y a random sample from a Poisson distribution.
#' @param shape the shape parameter of the \eqn{gamma} prior.
#' @param rate the rate parameter of the \eqn{gamma} prior. Note that the scale
#' is \eqn{1 / rate}
#' @param scale the scale parameter of the \eqn{gamma} prior
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced.
#' @param suppressOutput if \code{TRUE} then none of the output is printed to
#' console
#' @return An object of class 'Bolstad' is returned. This is a list with the
#' following components:
#' 
#' \item{prior}{the prior density assigned to \eqn{\mu}{mu}}
#' \item{likelihood}{the scaled likelihood function for \eqn{\mu}{mu} given
#' \eqn{y}} \item{posterior}{the posterior probability of \eqn{\mu}{mu} given
#' \eqn{y}} \item{shape}{the shape parameter for the \eqn{gamma} posterior}
#' \item{rate}{the rate parameter for the \eqn{gamma} posterior}
#' @seealso \code{\link{poisdp}} \code{\link{poisgcp}}
#' @keywords misc
#' @examples
#' 
#' ## simplest call with an observation of 4 and a gamma(1,1), i.e. an exponential prior on the
#' ## mu
#' poisgamp(4,1,1)
#' 
#' ##  Same as the previous example but a gamma(10,1) prior
#' poisgamp(4,10,1)
#' 
#' ##  Same as the previous example but an improper gamma(1,0) prior
#' poisgamp(4,1,0)
#' 
#' ## A random sample of 50 observations from a Poisson distribution with
#' ## parameter mu = 3 and  gamma(6,3) prior
#' y = rpois(50,3)
#' poisgamp(y,6,3)
#' 
#' ## In this example we have a random sample from a Poisson distribution
#' ## with an unknown mean. We will use a gamma(6,3) prior to obtain the
#' ## posterior gamma distribution, and use the R function qgamma to get a
#' ## 95% credible interval for mu
#' y = c(3,4,4,3,3,4,2,3,1,7)
#' results = poisgamp(y,6,3)
#' ci = qgamma(c(0.025,0.975),results$shape, results$rate)
#' cat(paste("95% credible interval for mu: [",round(ci[1],3), ",", round(ci[2],3)),"]\n")
#' 
#' ## In this example we have a random sample from a Poisson distribution
#' ## with an unknown mean. We will use a gamma(6,3) prior to obtain the
#' ## posterior gamma distribution, and use the R function qgamma to get a
#' ## 95% credible interval for mu
#' y = c(3,4,4,3,3,4,2,3,1,7)
#' results = poisgamp(y, 6, 3)
#' ci = quantile(results, c(0.025, 0.975))
#' cat(paste("95% credible interval for mu: [",round(ci[1],3), ",", round(ci[2],3)),"]\n")
#' 
#' 
#' @export poisgamp
poisgamp = function(y, shape, rate = 1, scale = 1 / rate,
                    plot = TRUE, suppressOutput = FALSE){
  n = length(y)
  y.sum = sum(y)

  if(is.null(y) || length(y)==0)
    stop("Error: y has no data")

  if(any(y < 0))
    stop("Error: y contains negative values")

  if(scale !=1 & rate == 1){
    rate = 1 / scale
  }
  
  r = shape
  v = rate
  
  if(r < 0 || v < 0)
    stop("Shape parameter and rate parameter must be greater than or equal to zero")

  if(!suppressOutput){
    cat("Summary statistics for data\n")
    cat("---------------------------\n")
    cat(paste("Number of observations:\t", n,"\n"))
    cat(paste("Sum of observations:\t", y.sum,"\n\n"))
  }
  
  if(v>0){                              ##proper gamma prior
    v.inv = 1/v
    k1 = qgamma(0.9999,r,v)
    k2 = k1/1000
    mu = seq(0, k1, by = k2)
    r1 = r + y.sum
    v1 = v + n
    v1.inv = 1 / v1

    prior = dgamma(mu,r, v)
    likelihood = matrix(0, ncol = length(mu), nrow = length(y))
    for(i in 1:length(mu)){
        likelihood[,i] = dpois(y,mu[i])
    }
    likelihood = apply(likelihood, 2, prod)
    posterior = dgamma(mu, r1, v1)

    k3 = qgamma(c(0.005,0.995), r1, v1)

    if(!suppressOutput){
      cat("Summary statistics for posterior\n")
      cat("--------------------------------\n")
      cat(paste("Shape parameter r:\t", r1,"\n"))
      cat(paste("Rate parameter v:\t",v1,"\n"))
      cat(paste("99% credible interval for mu:\t[",round(k3[1],2), ",",round(k3[2],2), "]\n"))
    }
    
    if(plot){
      y.max = max(prior,posterior)
      plot(mu,prior,ylim = c(0,1.1*y.max),xlab = expression(mu)
           ,ylab = "Density",
           ,main = "Shape of gamma prior and posterior\n for Poisson mean"
           ,type = "l",lty = 2,col = "red")
      lines(mu,posterior,lty = 3,col = "blue")
      legend("topleft", bty = "n", lty = 2:3, col = c("red","blue"),
             legend = c("Prior","Posterior"), cex = 0.7)
    }
  }else if(v == 0){
    r1 = r+y.sum
    v1 = v+n
    v1.inv = 1/v1
    k3 = qgamma(c(0.005,0.995),r1,v1)
    k4 = k3[2]/1000
    mu = seq(0,k3[2],by = k4)

    mu[1] = mu[2] ## fixes infinite upper bound problem

    if(!suppressOutput){
      cat("Summary statistics for posterior\n")
      cat("--------------------------------\n")
      cat(paste("Shape parameter r:\t", r1,"\n"))
      cat(paste("Rate parameter v:\t",v1,"\n"))
      cat(paste("99% credible interval :\t[",round(k3[1],2),", ",round(k3[2],2), "]\n"))
    }
    
    prior = mu^(r-1)
    kint = (2*sum(prior)-prior[1001])*k4/2
    prior = prior/kint

    likelihood = matrix(0, ncol = length(mu), nrow = length(y))
    for(i in 1:length(mu)){
        likelihood[,i] = dpois(y,mu[i])
    }
    likelihood = apply(likelihood, 2, prod)

    posterior = dgamma(mu, r1, v1)
    
    if(plot){
      y.max = max(prior,posterior)
      plot(mu,prior,ylim = c(0,1.1*y.max),xlab = expression(mu)
           ,ylab = "Density",
           ,main = "Shape of gamma prior and posterior\n for Poisson mean"
           ,type = "l",lty = 2,col = "red")
      lines(mu,posterior,lty = 3,col = "blue")
      legend("topleft", bty = "n", cex = 0.7,
             lty = 2:3, col = c("red", "blue"), legend = c("Prior","Posterior"))
    }
  }else{
    stop("Error: rate must be greater or equal to zero")
  }

  results = list(name = 'mu', param.x = mu, prior = prior, 
                 likelihood = likelihood, posterior = posterior,
                 mean = r1 / v1, 
                 var = r1 / v1^2,
                 cdf = function(m, ...){pgamma(m, shape = r1, rate = v1, ...)},
                 quantileFun = function(probs, ...){qgamma(probs, shape = r1, rate = v1, ...)},
                 mu = mu, # for backwards compatibility only
                 shape = r1, rate = v1)
  
  class(results) = 'Bolstad'
  invisible(results)
}



