#' Binomial sampling with a beta mixture prior
#' 
#' Evaluates and plots the posterior density for \eqn{\pi}{pi}, the probability
#' of a success in a Bernoulli trial, with binomial sampling when the prior
#' density for \eqn{\pi}{pi} is a mixture of two beta distributions,
#' \eqn{beta(a_0,b_0)} and \eqn{beta(a_1,b_1)}.
#' 
#' 
#' @param x the number of observed successes in the binomial experiment.
#' @param n the number of trials in the binomial experiment.
#' @param alpha0 a vector of length two containing the parameters,
#' \eqn{a_0}{a0} and \eqn{b_0}{b0}, for the first component beta prior - must
#' be greater than zero. By default the elements of alpha0 are set to 1.
#' @param alpha1 a vector of length two containing the parameters,
#' \eqn{a_1}{a1} and \eqn{b_1}{b1}, for the second component beta prior - must
#' be greater than zero. By default the elements of alpha1 are set to 1.
#' @param p The prior mixing proportion for the two component beta priors. That
#' is the prior is \eqn{p\times beta(a_0,b_0)+(1-p)\times
#' beta(a_1,b_1)}{p*beta(a0,b0)+(1-p)*beta(a1,b1)}. \eqn{p} is set to 0.5 by
#' default
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced
#' @return A list will be returned with the following components: \item{pi}{the
#' values of \eqn{\pi}{pi} for which the posterior density was evaluated}
#' \item{posterior}{the posterior density of \eqn{\pi}{pi} given \eqn{n} and
#' \eqn{x}} \item{likelihood}{the likelihood function for \eqn{\pi}{pi} given
#' \eqn{x} and \eqn{n}, i.e. the \eqn{binomial(n,\pi)}{binomial(n,pi)} density}
#' \item{prior}{the prior density of \eqn{\pi}{pi} density}
#' @seealso \code{\link{binodp}} \code{\link{binogcp}} \code{\link{normmixp}}
#' @keywords misc
#' @examples
#' 
#' ## simplest call with 6 successes observed in 8 trials and a 50:50 mix
#' ## of two beta(1,1) uniform priors
#' binomixp(6,8)
#' 
#' ## 6 successes observed in 8 trials and a 20:80 mix of a non-uniform
#' ## beta(0.5,6) prior and a uniform beta(1,1) prior
#' binomixp(6,8,alpha0=c(0.5,6),alpha1=c(1,1),p=0.2)
#' 
#' ## 4 successes observed in 12 trials with a 90:10 non uniform beta(3,3) prior
#' ## and a non uniform beta(4,12).
#' ## Plot the stored prior, likelihood and posterior
#' results = binomixp(4,12,c(3,3),c(4,12),0.9)$mix
#' 
#' par(mfrow=c(3,1))
#' y.lims = c(0,1.1*max(results$posterior,results$prior))
#' 
#' plot(results$pi,results$prior,ylim=y.lims,type="l"
#' 	,xlab=expression(pi),ylab="Density",main="Prior")
#' polygon(results$pi,results$prior,col="red")
#' 
#' plot(results$pi,results$likelihood,type="l"
#' 	,xlab=expression(pi),ylab="Density",main="Likelihood")
#' polygon(results$pi,results$likelihood,col="green")
#' 
#' plot(results$pi,results$posterior,ylim=y.lims,type="l"
#' 	,xlab=expression(pi),ylab="Density",main="Posterior")
#' polygon(results$pi,results$posterior,col="blue")
#' 
#' 
#' 
#' 
#' @export binomixp
binomixp = function(x, n, alpha0 = c(1,1), alpha1 = c(1,1), p = 0.5, plot = TRUE){

  if(n < x)
    stop("Error: n must be greater than or equal to x")

  if(length(alpha0) != 2 || length(alpha1) != 2)
    stop("Error: the parameters for the beta priors, alpha0 and alpha1, must have two elements each")

  if(any(alpha0 < 0) || any(alpha1 < 0))
    stop("Error: the parameters for the beta priors, alpha0 and alpha1, must be greater than zero")

  if(p <= 0 || p >= 1)
    stop("Error: the mixing proportion p must be in the interval (0,1) exclusive")


  i = 1:x
  log0 = sum(log(alpha0[1]+i-1)-log(alpha0[1]+alpha0[2]+i-1)+log(n-i+1)-log(i))
  log1 = sum(log(alpha1[1]+i-1)-log(alpha1[1]+alpha1[2]+i-1)+log(n-i+1)-log(i))
  i = (x+1):n
  log0 = log0+sum(log(alpha0[2]+i-x-1)-log(alpha0[1]+alpha0[2]+i-1))
  log1 = log1+sum(log(alpha1[2]+i-x-1)-log(alpha1[1]+alpha1[2]+i-1))

  f0 = exp(log0)
  f1 = exp(log1)

  cat("Prior probability of the data under component 0\n")
  cat("----------------------------\n")
  cat(paste("Log prob.:\t",signif(log0,3),"\nProbability:\t ",signif(f0,5),"\n\n"))

  cat("Prior probability of the data under component 1\n")
  cat("----------------------------\n")
  cat(paste("Log prob.:\t",signif(log1,3),"\nProbability:\t ",signif(f1,5),"\n\n"))


  q0 = p
  q1 = 1 - q0
  qp0 = q0 * f0 / (q0 * f0 + q1 * f1)
  qp1 = 1-qp0

  cat(paste("Post. mixing proportion for component 0:\t",signif(qp0, 3),"\n"))
  cat(paste("Post. mixing proportion for component 1:\t",signif(qp1, 3),"\n"))

  pi = seq(0.001, 0.999, by = 0.001)

  prior.0 = dbeta(pi, alpha0[1], alpha0[2])
  prior.1 = dbeta(pi, alpha1[1], alpha1[2])
  prior = q0 * prior.0 + q1 * prior.1

  alpha0.post = alpha0 + c(x, n - x)
  alpha1.post = alpha1 + c(x, n - x)

  posterior.0 = dbeta(pi, alpha0.post[1], alpha0.post[2])
  posterior.1 = dbeta(pi, alpha1.post[1], alpha1.post[2])
  posterior = qp0 * posterior.0 + qp1 * posterior.1

  loglik = x * log(pi) + (n - x) * log(1 - pi)
  loglik = loglik - max(loglik)
  likelihood = exp(loglik)

  normalizing.factor = sum(likelihood) / length(likelihood)
  likelihood = likelihood/normalizing.factor

  if(plot){
    o.par = par(mfrow=c(2, 2))
  
    ##plot the priors and the mixture prior
  
    y.max = max(prior.0,prior.1,prior)
  
    plot(pi,prior.0,ylim=c(0,y.max*1.1),
         xlab=expression(pi),ylab="Density"
         ,main="Mixture prior and its components"
         ,type="l",lty=2,col="black")
    lines(pi,prior.1,lty=3,col="red")
    lines(pi,prior,lty=1,col="green")
    legend("topleft", cex = 0.7, bty = "n", 
           legend=c(expression(prior[0]),expression(prior[1])
                     ,expression(prior[mix])),lty=c(2,3,1),col=c("black","red","green"))
  
    ##plot the posteriors and the mixture posterior
  
    y.max = max(posterior.0,posterior.1,posterior)
  
    plot(pi,posterior.0,ylim=c(0,y.max*1.1),
         xlab=expression(pi),ylab="Density"
         ,main="Mixture posterior and its components"
         ,type="l",lty=2,col="black")
    lines(pi,posterior.1,lty=3,col="red")
    lines(pi,posterior,lty=1,col="green")
    legend("topleft", bty = "n",
           legend=c(expression(posterior[0]),expression(posterior[1])
                     ,expression(posterior[mix])),
           lty = c(2, 3, 1), col=c("black", "red", "green"))
  
    ##plot the mixture posterior likelihood and mixture posterior
  
    y.max = max(prior,posterior,likelihood)
  
    plot(pi,prior,ylim=c(0,y.max*1.1),
         xlab = expression(pi),ylab="Density"
         ,main = "Mixture prior, likelihood and mixture posterior"
         ,type = "l",lty=2,col="black")
    lines(pi,likelihood,lty=3,col="red")
    lines(pi,posterior,lty=1,col="green")
    legend("topleft", bty = "n",
           legend=c(expression(prior[mix]),expression(likelihood)
                    ,expression(posterior[mix])),
           lty=c(2,3,1), col = c("black", "red", "green"))
  
    par(o.par)
  }
  
  results.comp1 = list(name = 'pi', param.x = pi, prior = prior.0,
                    likelihood = likelihood, posterior = posterior.0)
  class(results.comp1) = "Bolstad"
  
  results.comp2 = list(name = 'pi', param.x = pi, prior = prior.1,
                    likelihood = likelihood, posterior = posterior.1)
  class(results.comp2) = "Bolstad"
  
  results.mix = list(name = 'pi', param.x = pi ,prior = prior,
                 likelihood = likelihood, posterior = posterior,
                 pi = pi #for backwards compat. only
                 )
  class(results.mix) = 'Bolstad'
  invisible(list(comp1 = results.comp1, comp2 = results.comp2, mix = results.mix))
}





