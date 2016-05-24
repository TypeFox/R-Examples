#' Bayesian inference on a normal mean with a discrete prior
#' 
#' Evaluates and plots the posterior density for \eqn{\mu}{mu}, the mean of a
#' normal distribution, with a discrete prior on \eqn{\mu}{mu}
#' 
#' 
#' @param x a vector of observations from a normal distribution with unknown
#' mean and known std. deviation.
#' @param sigma.x the population std. deviation of the normal distribution
#' @param mu a vector of possibilities for the probability of success in a
#' single trial. If mu is NULL then a uniform prior is used.
#' @param mu.prior the associated prior probability mass.
#' @param n.mu the number of possible \eqn{\mu}{mu} values in the prior
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced.
#' @return A list will be returned with the following components: \item{mu}{the
#' vector of possible \eqn{\mu}{mu} values used in the prior}
#' \item{mu.prior}{the associated probability mass for the values in
#' \eqn{\mu}{mu}} \item{likelihood}{the scaled likelihood function for
#' \eqn{\mu}{mu} given \eqn{x} and \eqn{\sigma_x}{sigma.x}}
#' \item{posterior}{the posterior probability of \eqn{\mu}{mu} given \eqn{x}
#' and \eqn{\sigma_x}{sigma.x}}
#' @seealso \code{\link{normnp}} \code{\link{normgcp}}
#' @keywords misc
#' @examples
#' 
#' ## generate a sample of 20 observations from a N(-0.5,1) population
#' x = rnorm(20,-0.5,1)
#' 
#' ## find the posterior density with a uniform prior on mu
#' normdp(x,1)
#' 
#' ## find the posterior density with a non-uniform prior on mu
#' mu = seq(-3,3,by=0.1)
#' mu.prior = runif(length(mu))
#' mu.prior = sort(mu.prior/sum(mu.prior))
#' normdp(x,1,mu,mu.prior)
#' 
#' ## Let mu have the discrete distribution with 5 possible
#' ## values, 2, 2.5, 3, 3.5 and 4, and associated prior probability of
#' ## 0.1, 0.2, 0.4, 0.2, 0.1 respectively. Find the posterior
#' ## distribution after a drawing random sample of n = 5 observations
#' ## from a N(mu,1) distribution y = [1.52, 0.02, 3.35, 3.49, 1.82]
#' mu = seq(2,4,by=0.5)
#' mu.prior = c(0.1,0.2,0.4,0.2,0.1)
#' y = c(1.52,0.02,3.35,3.49,1.82)
#' normdp(y,1,mu,mu.prior)
#' 
#' @export normdp
normdp = function(x, sigma.x = NULL, mu = NULL, mu.prior = NULL, n.mu = 50, plot = TRUE){

  ## x - the vector of observations
  ## sigma.x - the population standard deviation
  ## mu - vector of possible values of the population mean
  ## mu.prior - the associated prior probability mass
  ## n.mu - if mu is NULL then a uniform prior with n.mu points is used
  ## the likelihood and posterior are returned as a
  ## list

  if(n.mu < 3)
    stop("Number of prior values of theta must be greater than 2")

  if(is.null(mu)){
    mu = seq(min(x)-sigma.x,max(x)+sigma.x,length = n.mu)
    mu.prior = rep(1/n.mu,n.mu)
  }

  mx = mean(x)

  if(is.null(sigma.x)){
    sigma.x = sd(x-mx)
    cat(paste("Standard deviation of the residuals :"
              ,signif(sigma.x,4),"\n",sep=""))
  }else{
    if(sigma.x>0){
      cat(paste("Known standard deviation :",signif(sigma.x,4),"\n",sep=""))
    }else{
      stop("The standard deviation must be greater than zero")
    }
  }

  if(any(mu.prior<0) | any(mu.prior>1))
    stop("Prior probabilities must be between 0 and 1 inclusive")

  if(round(sum(mu.prior),7)!=1){
    warning("The prior probabilities did not sum to 1, therefore the prior has been normalized")
    mu.prior = mu.prior/sum(mu.prior)
  }

  n.mu = length(mu)
  nx = length(x)
  snx = sigma.x^2/nx
  likelihood = exp(-0.5*(mx-mu)^2/snx)
  posterior = likelihood*mu.prior/sum(likelihood*mu.prior)

  if(plot){
    plot(mu,posterior,ylim=c(0,1.1*max(posterior,mu.prior)),pch = 20, col = "blue",
         xlab=expression(mu),ylab=expression(Probabilty(mu)))
    points(mu, mu.prior, pch = 20, col = "red")
  
    legend("topleft", bty = "n", fill = c("blue", "red"),
           legend = c("Posterior","Prior"), cex = 0.7)
  }
  mx = sum(mu * posterior)
  vx = sum((mu - mx)^2 * posterior)
  
  results = list(name = 'mu', param.x = mu.prior, 
                 prior = mu.prior, 
                 likelihood = likelihood, posterior = posterior,
                 mean = mx, var = vx,
                 cdf = function(x, ...)cumDistFun(x, mu, posterior),
                 quantileFun = function(probs, ...)qFun(probs, mu, posterior))
  class(results) = 'Bolstad'
  
  invisible(results)
}
