#' Binomial sampling with a beta prior
#' 
#' Evaluates and plots the posterior density for \eqn{\pi}{pi}, the probability
#' of a success in a Bernoulli trial, with binomial sampling and a continous
#' \eqn{beta(a,b)} prior.
#' 
#' 
#' @param x the number of observed successes in the binomial experiment.
#' @param n the number of trials in the binomial experiment.
#' @param a parameter for the beta prior - must be greater than zero
#' @param b parameter for the beta prior - must be greater than zero
#' @param pi A rannge of values for the prior to be calculated over.
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced.
#' @return An object of class 'Bolstad' is returned. This is a list with the
#' following components: \item{prior}{the prior density of \eqn{\pi}{pi}, i.e.
#' the \eqn{beta(a,b)} density} \item{likelihood}{the likelihood of \eqn{x}
#' given \eqn{\pi}{pi} and \eqn{n}, i.e. the
#' \eqn{binomial(n,\pi)}{binomial(n,pi)} density} \item{posterior}{the
#' posterior density of \eqn{\pi}{pi} given \eqn{x} and \eqn{n} - i.e. the
#' \eqn{beta(a+x,b+n-x)} density} \item{pi}{the values of \eqn{\pi}{pi} for
#' which the posterior density was evaluated} \item{mean}{the posterior mean}
#' \item{var}{the posterior variance} \item{sd}{the posterior std. deviation}
#' \item{quantiles}{a set of quantiles from the posterior} \item{cdf}{a
#' cumulative distribution function for the posterior} \item{quantileFun}{a
#' quantile function for the posterior}
#' @seealso \code{\link{binodp}} \code{\link{binogcp}}
#' @keywords misc
#' @examples
#' 
#' ## simplest call with 6 successes observed in 8 trials and a beta(1,1) uniform
#' ## prior
#' binobp(6,8)
#' 
#' ## 6 successes observed in 8 trials and a non-uniform beta(0.5,6) prior
#' binobp(6,8,0.5,6)
#' 
#' ## 4 successes observed in 12 trials with a non uniform beta(3,3) prior
#' ## plot the stored prior, likelihood and posterior
#' results = binobp(4, 12, 3, 3)
#' decomp(results)
#' 
#' 
#' @export binobp
binobp = function(x, n, a = 1, b = 1, pi = seq(0.01, 0.999, by = 0.001), plot = TRUE){
  
  ## n - the number of trials in the binomial
  ## x - the number of observed successes
  ## a,b  - the parameters of the Beta prior density (must be > 0)
  ## the prior, likelihood, posterior, mean, variance and
  ## std. deviation are returned as a list
  
  if(x > n)
    stop("The number of observed successes (x) must be smaller than the number of trials (n)")
  
  if(a <= 0 || b <= 0)
    stop("The parameters of the prior must be greater than zero")
  
  
  prior = dbeta(pi, a, b)
  likelihood = dbinom(x, n, prob = pi)
  posterior = dbeta(pi, a + x, b + n - x)
  
  if(plot){
    plot(posterior ~ pi, ylim = c(0, 1.1 * max(posterior, prior)), type="l",
         lty = 1,
         xlab = expression(pi),
         ylab = "Density",
         col = "blue")
    lines(prior ~ pi, lty = 2, col = "red")
    ## left = min(pi) + diff(range(pi)) * 0.05
    legend("topleft", bty = "n", lty = 1:2, legend = c("Posterior","Prior"),
           col = c("blue","red"), cex = 0.7)
  }
  
  m1 = (a + x) / (a + b + n)
  v1 = m1 * (1 - m1) / (a + b + n + 1)
  s1 = sqrt(v1)
  
  cat(paste("Posterior Mean           : ",round(m1,7),"\n"))
  cat(paste("Posterior Variance       : ",round(v1,7),"\n"))
  cat(paste("Posterior Std. Deviation : ",round(s1,7),"\n"))
  
  probs = c(0.005,0.01,0.025,0.05,0.5,0.95,0.975,0.99,0.995)
  qtls = qbeta(probs,a+x,b+n-x)
  names(qtls) = probs
  
  cat("\nProb.\tQuantile \n")
  cat("------\t---------\n")
  for(i in 1:length(probs)){
    cat(sprintf("%5.3f\t%9.7f\n", round(probs[i],3),round(qtls[i],7)))
  }
  
  
  results = list(name = 'pi', param.x = pi, prior = prior, likelihood = likelihood, posterior = posterior,
                 pi = pi, # for backwards compat. only
                 mean = m1, var = v1, sd = s1, quantiles = qtls,
                 cdf = function(y,...){pbeta(y, shape1 = a + x, shape2 = b + n - x, ...)},
                 quantileFun = function(probs, ...){
                   qbeta(probs, shape1 = a + x, shape2 = b + n - x, ...)})
  
  class(results) = 'Bolstad'
  invisible(results)
}
