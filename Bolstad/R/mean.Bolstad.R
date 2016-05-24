#'Calculate the posterior mean
#'
#'Calculate the posterior mean of an object of class \code{Bolstad}. If the 
#'object has a member \code{mean} then it will return this value otherwise it 
#'will calculate \eqn{\int_{-\infty}^{+\infty}\theta f(\theta|x).d\theta} using 
#'linear interpolation to approximate the density function and numerical 
#'integration where \eqn{\theta} is the variable for which we want to do 
#'Bayesian inference, and \eqn{x} is the data.
#'
#'@param x An object of class \code{Bolstad}
#'@param \dots Any other arguments. This parameter is currently ignored but it 
#'  could be useful in the future to deal with problematic data.
#'@return The posterior mean of the variable of inference given the data.
#'@examples 
#' # The useful of this method is really highlighted when we have a general 
#' # continuous prior. In this example we are interested in the posterior mean of 
#' # an normal mean. Our prior is triangular over [-3, 3]
#'set.seed(123)
#'x = rnorm(20, -0.5, 1)
#'mu = seq(-3, 3, by = 0.001)
#'mu.prior = rep(0, length(mu))
#'mu.prior[mu <= 0] = 1 / 3 + mu[mu <= 0] / 9
#'mu.prior[mu > 0] = 1 / 3 - mu[mu > 0] / 9
#'results = normgcp(x, 1, density = "user", mu = mu, mu.prior = mu.prior)
#'mean(results)
#'@method mean Bolstad
#'@export
mean.Bolstad = function(x, ...){
  if(any(grepl("mean", names(x))))
    return(x$mean)
  
  xVals = x$param.x
  fx = approxfun(xVals, xVals * x$posterior)
  
  return(integrate(fx, min(xVals), max(xVals))$value)
}

