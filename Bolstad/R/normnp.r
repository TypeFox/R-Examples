#' Bayesian inference on a normal mean with a normal prior
#' 
#' Evaluates and plots the posterior density for \eqn{\mu}{mu}, the mean of a
#' normal distribution, with a normal prior on \eqn{\mu}{mu}
#' 
#' 
#' @param x a vector of observations from a normal distribution with unknown
#' mean and known std. deviation.
#' @param m.x the mean of the normal prior
#' @param s.x the standard deviation of the normal prior
#' @param sigma.x the population std. deviation of the normal distribution. If
#' this value is NULL, which it is by default, then a flat prior is used and
#' m.x and s.x are ignored
#' @param mu a vector of prior possibilities for the true mean. If this is \code{null},
#' then a set of values centered on the sample mean is used.
#' @param n.mu the number of possible \eqn{\mu}{mu} values in the prior
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced.
#' @return A list will be returned with the following components: \item{mu}{the
#' vector of possible \eqn{\mu}{mu} values used in the prior}
#' \item{mu.prior}{the associated probability mass for the values in
#' \eqn{\mu}{mu}} \item{likelihood}{the scaled likelihood function for
#' \eqn{\mu}{mu} given \eqn{x} and \eqn{\sigma_x}{sigma.x}}
#' \item{posterior}{the posterior probability of \eqn{\mu}{mu} given \eqn{x}
#' and \eqn{\sigma_x}{sigma.x}} \item{mean}{the posterior mean} \item{sd}{the
#' posterior standard deviation} \item{qtls}{a selection of quantiles from the
#' posterior density}
#' @seealso \code{\link{normdp}} \code{\link{normgcp}}
#' @keywords misc
#' @examples
#' 
#' ## generate a sample of 20 observations from a N(-0.5,1) population
#' x = rnorm(20,-0.5,1)
#' 
#' ## find the posterior density with a N(0,1) prior on mu
#' normnp(x,sigma=1)
#' 
#' ## find the posterior density with N(0.5,3) prior on mu
#' normnp(x,0.5,3,1)
#' 
#' ## Find the posterior density for mu, given a random sample of 4
#' ## observations from N(mu,sigma^2=1), y = [2.99, 5.56, 2.83, 3.47],
#' ## and a N(3,sd=2)$ prior for mu
#' y = c(2.99,5.56,2.83,3.47)
#' normnp(y,3,2,1)
#' 
#' @export normnp
normnp = function(x, m.x = 0 , s.x = 1, sigma.x = NULL, mu = NULL, n.mu = max(100, length(mu)), 
                  plot = TRUE){

  mean.x = mean(x)
  n.x = length(x)

  if(is.null(sigma.x)){
    sigma.x = sd(x - mean.x)
    cat(paste("Standard deviation of the residuals :",signif(sigma.x,4),"\n",sep=""))
  }else{
    if(sigma.x > 0){
      cat(paste("Known standard deviation :",signif(sigma.x,4),"\n",sep=""))
    }else{
      stop("Standard deviation sigma.x must be greate than zero")
    }
  }

  if(is.null(mu)){
    lb = ub = 0
    if(s.x <= 0){
      lb = mean.x - 3.5 * sigma.x / sqrt(n.x)
      ub = mean.x + 3.5 * sigma.x / sqrt(n.x)
    }else{
      lb = m.x - 3.5 * s.x
      ub = m.x + 3.5 * s.x
    }
    mu = seq(lb, ub, length = n.mu)
  }else{
    if(length(mu) < n.mu)
      mu = seq(min(mu), max(mu), length = n.mu)
  }
  
  if(s.x <= 0){
    prior.precision = 0
    m.x = 0
    lb = min(mu)
    ub = max(mu)  
    mu.prior = rep(1 / (ub - lb), n.mu)
  }else{
    mu.prior = dnorm(mu, m.x, s.x)
    prior.precision = 1/s.x^2
  }


  likelihood = exp(-n.x/(2*sigma.x^2)*(mean.x-mu)^2)

  post.precision = prior.precision+(n.x/sigma.x^2)
  post.sd = sqrt(1/post.precision)
  post.mean = (prior.precision/post.precision*m.x)+((n.x/sigma.x^2)/post.precision*mean.x)

  cat(paste("Posterior mean           : ",round(post.mean,7),"\n",sep=""))
  cat(paste("Posterior std. deviation : ",round(post.sd,7),"\n",sep=""))

  posterior = dnorm(mu,post.mean,post.sd)

  if(plot){
    plot(mu,posterior,ylim=c(0,1.1*max(posterior,mu.prior)),type="l",
         lty=1,col="blue",
         xlab=expression(mu),ylab=expression(Probabilty(mu)),
         main="Shape of prior and posterior")
    lines(mu,mu.prior,lty=2,col="red")
  
    left = min(mu)+diff(range(mu))*0.05
    legend("topleft",bty = "n", lty = 1:2, col = c("blue", "red"),
           legend = c("Posterior", "Prior"), cex = 0.7)
  }
  probs = c(0.005,0.01,0.025,0.05,0.5,0.95,0.975,0.99,0.995)
  qtls = qnorm(probs,post.mean,post.sd)
  names(qtls) = probs

  cat("\nProb.\tQuantile \n")
  cat("------\t----------\n")
  for(i in 1:length(probs)){
    cat(sprintf("%5.3f\t%10.7f\n", round(probs[i],3), round(qtls[i],7)))
  }

  results = list(name = 'mu', param.x = mu, prior = mu.prior, likelihood = likelihood, posterior = posterior,
                 mean = post.mean, 
                 var = post.sd^2, 
                 sd = post.sd, quantiles = qtls,
                 cdf = function(x, ...)pnorm(x, post.mean, post.sd, ...),
                 quantileFun = function(probs, ...)qnorm(probs, post.mean, post.sd, ...))
  class(results) = 'Bolstad'
  
  invisible(results)
}
