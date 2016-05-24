#' Bayesian inference on a normal mean with a mixture of normal priors
#' 
#' Evaluates and plots the posterior density for \eqn{\mu}{mu}, the mean of a
#' normal distribution, with a mixture of normal priors on \eqn{\mu}{mu}
#' 
#' 
#' @param x a vector of observations from a normal distribution with unknown
#' mean and known std. deviation.
#' @param sigma.x the population std. deviation of the observations.
#' @param prior0 the vector of length 2 which contains the means and standard
#' deviation of your precise prior.
#' @param prior1 the vector of length 2 which contains the means and standard
#' deviation of your vague prior.
#' @param mu a vector of prior possibilities for the mean. If it is \code{NULL},
#' then a vector centered on the sample mean is created.
#' @param n.mu the number of possible \eqn{\mu}{mu} values in the prior.
#' @param p the mixing proportion for the two component normal priors.
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced.
#' @return A list will be returned with the following components: \item{mu}{the
#' vector of possible \eqn{\mu}{mu} values used in the prior} \item{prior}{the
#' associated probability mass for the values in \eqn{\mu}{mu}}
#' \item{likelihood}{the scaled likelihood function for \eqn{\mu}{mu} given
#' \eqn{x} and \eqn{\sigma_x}{sigma.x}} \item{posterior}{the posterior
#' probability of \eqn{\mu}{mu} given \eqn{x} and \eqn{\sigma_x}{sigma.x}}
#' @seealso \code{\link{binomixp}} \code{\link{normdp}} \code{\link{normgcp}}
#' @keywords misc
#' @examples
#' 
#' ## generate a sample of 20 observations from a N(-0.5, 1) population
#' x = rnorm(20, -0.5, 1)
#' 
#' ## find the posterior density with a N(0, 1) prior on mu - a 50:50 mix of
#' ## two N(0, 1) densities
#' normmixp(x, 1, c(0, 1), c(0, 1))
#' 
#' ## find the posterior density with 50:50 mix of a N(0.5, 3) prior and a
#' ## N(0, 1) prior on mu
#' normmixp(x, 1, c(0.5, 3), c(0, 1))
#' 
#' ## Find the posterior density for mu, given a random sample of 4
#' ## observations from N(mu, 1), y = [2.99, 5.56, 2.83, 3.47],
#' ## and a 80:20 mix of a N(3, 2) prior and a N(0, 100) prior for mu
#' x = c(2.99, 5.56, 2.83, 3.47)
#' normmixp(x, 1, c(3, 2), c(0, 100), 0.8)
#' 
#' @export normmixp
normmixp = function(x, sigma.x, prior0, prior1, p = 0.5, mu = NULL, n.mu = max(100, length(mu)),
                    plot = TRUE){

  if(length(x) == 0)
    stop("Error: x must contain at least one observation")

  if(sigma.x <= 0)
    stop("Error: sigma.x must be greater than zero")

  if(length(prior0) != 2 || length(prior1) != 2)
    stop("Error: there must be 2 parameters for each prior, 2 means and 2 standard deviations")

  prior.means = c(prior0[1], prior1[1])
  prior.sds = c(prior0[2], prior1[2])

  if(any(prior.sds <= 0))
    stop("Error: the prior standard deviations should be greater than zero")

  if(p <= 0 || p >= 1)
    stop("Error: p should be between 0 and 1 exclusive")

  n = length(x)
  x.bar = mean(x)
  q0 = p
  q1 = 1 - p

  post.prec0 = (1 / prior.sds[1]^2) + (n / sigma.x^2)
  post.var0 = 1 / post.prec0
  post.sd0 = sqrt(post.var0)
  post.mean0 = (prior.means[1] / (prior.sds[1]^2 * post.prec0)) + (n * x.bar / (sigma.x^2 * post.prec0))


  post.prec1 = (1 / prior.sds[2]^2) + (n / sigma.x^2)
  post.var1 = 1 / post.prec1
  post.sd1 = sqrt(post.var1)
  post.mean1 = (prior.means[2] / (prior.sds[2]^2 * post.prec1)) + (n * x.bar / (sigma.x^2 * post.prec1))

  cat("Posterior summary statistics of component 0\n")
  cat("--------------------------------------------\n")
  cat(paste("Mean:\t\t", signif(post.mean0, 3), "\n"))
  cat(paste("Std. Dev.:\t", signif(post.sd0, 4), "\n"))
  cat(paste("Variance:\t", signif(post.var0, 4), "\n"))
  cat(paste("Precision:\t", signif(post.prec0, 4), "\n\n"))

  cat("Posterior summary statistics of component 1\n")
  cat("-------------------------------------------\n")
  cat(paste("Mean:\t\t", signif(post.mean1, 3), "\n"))
  cat(paste("Std. Dev.:\t", signif(post.sd1, 4), "\n"))
  cat(paste("Variance:\t", signif(post.var1, 4), "\n"))
  cat(paste("Precision:\t", signif(post.prec1, 4), "\n\n"))


  sd.x = sqrt(sigma.x^2 / n + prior.sds[1]^2)
  f0 = dnorm(x.bar, prior.means[1], sd.x)

  cat("Predictive density of the sample mean under component 0\n")
  cat("------------------------------------------------------\n")
  cat(paste("Sample mean:\t", signif(x.bar, 3), "\n"))
  cat(paste("Pred. mean:\t", signif(prior.means[1], 3), "\n"))
  cat(paste("Pred. SD:\t", signif(sd.x, 4), "\n"))
  cat(paste("Density:\t", signif(f0, 4), "\n\n"))


  sd.x = sqrt(sigma.x^2 / n + prior.sds[2]^2)
  f1 = dnorm(x.bar, prior.means[2], sd.x)

  cat("Predictive density of the sample mean under component 1\n")
  cat("------------------------------------------------------\n")
  cat(paste("Sample mean:\t", signif(x.bar, 3), "\n"))
  cat(paste("Pred. mean:\t", signif(prior.means[2], 3), "\n"))
  cat(paste("Pred. SD:\t", signif(sd.x, 4), "\n"))
  cat(paste("Density:\t", signif(f1, 4), "\n\n"))

  qp0 = q0 * f0 / (q0 * f0 + q1 * f1)
  qp1 = 1 - qp0

  cat(paste("Post. mixing proportion for component 0:\t", signif(qp0, 3), "\n"))
  cat(paste("Post. mixing proportion for component 1:\t", signif(qp1, 3), "\n"))

  step.size = k1 = k2 = 0
  if(is.null(mu)){
    k1 = qnorm(0.0001, prior0[1], prior0[2])
    k2 = qnorm(0.9999, prior0[1], prior0[2])
    step.size = (k2 - k1) / 1000
    mu = seq(k1, k2, by = step.size)
  }
  
  if(length(mu) < n.mu){
    k1 = min(mu)
    k2 = max(mu)
    mu = seq(k1, k2, length = n.mu)
    step.size = diff(mu)[1]
  }else{
    k1 = min(mu)
    k2 = max(mu)
    step.size = diff(mu)[1]
  }
  
  
  prior.0 = dnorm(mu, prior.means[1], prior.sds[1])
  prior.1 = dnorm(mu, prior.means[2], prior.sds[2])
  prior = q0 * prior.0 + q1 * prior.1

  posterior.0 = dnorm(mu, post.mean0, post.sd0)
  posterior.1 = dnorm(mu, post.mean1, post.sd1)
  posterior = qp0 * posterior.0 + qp1 * posterior.1

  loglik = -(mu - x.bar)^2 / (2 * sigma.x^2 / n)
  loglik = loglik - max(loglik)
  likelihood = exp(loglik)

  normalizing.factor = sum(likelihood) * step.size
  likelihood = likelihood / normalizing.factor

  f.mu = approxfun(mu, likelihood)
  cat(paste("\nIntegral of likelihood over mu: ", round(integrate(f.mu, k1, k2)$value, 5), "\n"))

  if(plot){
    o.par = par(mfrow = c(2, 2))
  
    ##plot the priors and the mixture prior
  
    y.max = max(prior.0, prior.1, prior)
  
    plot(mu, prior.0, ylim = c(0, y.max * 1.1),
         xlab = expression(mu), ylab = "Density"
         , main = "Mixture prior and it's components"
         , type = "l", lty = 2, col = "black")
    lines(mu, prior.1, lty = 3, col = "red")
    lines(mu, prior, lty = 1, col = "blue")
    legend("topleft", bty = "n", cex = 0.7,
           legend = c(expression(prior[0]), expression(prior[1])
                      , expression(prior[mix])),
           lty = c(2, 3, 1), col = c("black", "red", "blue"))
  
    ##plot the posteriors and the mixture posterior
  
    y.max = max(posterior.0, posterior.1, posterior)
  
    plot(mu, posterior.0, ylim = c(0, y.max * 1.1),
         xlab = expression(mu), ylab = "Density"
         , main = "Mixture posterior and it's components"
         , type = "l", lty = 2, col = "black")
    lines(mu, posterior.1, lty = 3, col = "red")
    lines(mu, posterior, lty = 1, col = "blue")
    legend("topleft", bty = "n", cex = 0.7, 
           legend = c(expression(posterior[0]), expression(posterior[1])
                      , expression(posterior[mix])),
           lty = c(2, 3, 1), col = c("black", "red", "blue"))
    ##plot the mixture posterior likelihood and mixture posterior
  
    y.max = max(prior, posterior, likelihood)
  
    plot(mu, prior, ylim = c(0, y.max * 1.1),
         xlab = expression(mu), ylab = "Density"
         , main = "Mixture prior, likelihood and mixture posterior"
         , type = "l", lty = 2, col = "black")
    lines(mu, likelihood, lty = 3, col = "red")
    lines(mu, posterior, lty = 1, col = "blue")
    legend("topleft", bty = "n", cex = 0.7,
           legend = c(expression(prior[mix]), expression(likelihood)
                      , expression(posterior[mix])),
           lty = c(2, 3, 1), col = c("black", "red", "blue"))
  
    par(o.par)
  }
  comp1 = list(name = 'mu', param.x =  mu, 
               prior = prior.0, likelihood = likelihood, posterior = posterior.0,
               mean  = post.mean0, var = post.var0,
               cdf = function(q, ...){pnorm(q, mean = post.mean0, sd = post.sd0, ...)},
               quantileFun = function(probs, ...){qnorm(probs, mean = post.mean0, sd = post.sd0, ...)})
  class(comp1) = "Bolstad"
  
  comp2 = list(name = 'mu', param.x =  mu, 
               prior = prior.1, likelihood = likelihood, posterior = posterior.1,
               mean  = post.mean1, var = post.var1,
               cdf = function(q, ...){pnorm(q, mean = post.mean1, sd = post.sd1, ...)},
               quantileFun = function(probs, ...){qnorm(probs, mean = post.mean1, sd = post.sd1, ...)})
  class(comp2) = "Bolstad"
  
  mix = list(name = 'mu', param.x = mu, 
             prior = prior, likelihood = likelihood, posterior = posterior,
             mu = mu ## for backwards compatibility only
  )
  class(mix) = "Bolstad"
  invisible(list(comp1 = comp1, comp2 = comp2, mix = mix))
}





