#' Bayesian inference for a normal standard deviation with a scaled inverse
#' chi-squared distribution
#' 
#' Evaluates and plots the posterior density for \eqn{\sigma}{sigma}, the
#' standard deviation of a Normal distribution where the mean \eqn{\mu}{mu} is
#' known
#' 
#' 
#' @param y a random sample from a
#' \eqn{normal(\mu,\sigma^2)}{normal(mu,sigma^2)} distribution.
#' @param mu the known population mean of the random sample.
#' @param S0 the prior scaling factor.
#' @param kappa the degrees of freedom of the prior.
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced.
#' @param \dots this allows the arguments \code{cred.int} (which is logical), and 
#' \code{alpha} (numerical between 0 and 1 exclusive) to be specified for compatibility
#' with previous versions. A warning will be issued about these arguments being 
#' deprecated which is why there are no examples using them.
#' @return A list will be returned with the following components:
#' 
#' \item{sigma}{the vaules of \eqn{\sigma}{sigma} for which the prior,
#' likelihood and posterior have been calculated} \item{prior}{the prior
#' density for \eqn{\sigma}{sigma}} \item{likelihood}{the likelihood function
#' for \eqn{\sigma}{sigma} given \eqn{y}{y}} \item{posterior}{the posterior
#' density of \eqn{\mu}{sigma} given \eqn{y}{y}} \item{S1}{the posterior
#' scaling constant} \item{kappa1}{the posterior degrees of freedom}
#' @keywords misc
#' @examples
#' 
#' ## Suppose we have five observations from a normal(mu, sigma^2)
#' ## distribution mu = 200 which are 206.4, 197.4, 212.7, 208.5.
#' y = c(206.4, 197.4, 212.7, 208.5, 203.4)
#' 
#' ## We wish to choose a prior that has a median of 8. This happens when
#' ## S0 = 29.11 and kappa = 1
#' nvaricp(y,200,29.11,1)
#' 
#' ##  Same as the previous example but a calculate a 95% credible
#' ## interval for sigma. NOTE this method has changed
#' results = nvaricp(y,200,29.11,1)
#' quantile(results, probs = c(0.025, 0.975))
#' @export nvaricp
nvaricp = function(y, mu, S0, kappa, plot = TRUE, ...){
  
  dots = list(...)
  cred.int = dots[[pmatch("cred.int", names(dots))]]
  alpha = dots[[pmatch("alpha", names(dots))]]
  
  if(is.null(cred.int))
    cred.int = FALSE
  
  if(is.null(alpha))
    alpha = 0.05
  
  n = length(y)
  SST = sum((y-mu)^2)

  if(kappa > 0){
    S1 = S0 + SST
    kappa1 = kappa + n

    k1 = qchisq(0.01, kappa)
    k2 = S0/k1
    k3 = sqrt(k2)


    sigma = seq(0,k3,length = 1002)[-1]
    k4 = diff(sigma)[1]
    sigma.sq = sigma^2
    log.prior =  -((kappa-1)/2+1)*log(sigma.sq)-S0/(2*sigma.sq)
    prior = exp(log.prior)

    log.like =  -(n/2)*log(sigma.sq)-SST/(2*sigma.sq)
    likelihood =  exp(log.like)


    kint = ((2*sum(prior))-prior[1]-prior[1001])*k4/2*0.99
    prior = prior/kint

    posterior = prior*likelihood
    kint = ((2*sum(posterior))-posterior[1]-posterior[1001])*k4/2
    posterior = posterior/kint

    y.max = max(c(prior,posterior))
    k1 = qchisq(0.01,kappa1)
    k2 = S1/k1
    k3 = sqrt(k2)

    if(plot){
      plot(sigma, prior, type = "l", col = "blue", ylim = c(0, 1.1 * y.max), xlim = c(0,k3),
           main = expression(paste("Shape of Inverse ", chi^2," and posterior for ", sigma, sep = "")),
           xlab = expression(sigma),
           ylab = "Density")
      lines(sigma, posterior, lty = 1, col = "red")
      legend("topleft", lty = 1, lwd = 2, col = c("blue","red"), legend = c("Prior", "Posterior"), bty = "n")
    }
  }else if(kappa == 0){                   ## Jeffrey's prior
    S = 0
    S1 = S+SST
    kappa1 = kappa+n

    k1 = qchisq(0.001,kappa1)
    k2 = S1/k1
    k3 = sqrt(k2)
    k4 = k3/1000

    sigma = seq(0,k3,length = 1002)[-1]
    sigma.sq = sigma^2
    likelihood = NULL

    log.posterior =  -((kappa1-1)/2+1)*log(sigma.sq)-S1/(2*sigma.sq)
    posterior = exp(log.posterior)
    kint = ((2*sum(posterior))-posterior[1]-posterior[1001])*k4/(2*.999)
    posterior = posterior/kint

    log.prior =  -((kappa-1)/2+1)*log(sigma.sq)-S0/(2*sigma.sq)
    prior = exp(log.prior)
    kint = ((2*sum(prior))-prior[1]-prior[1001])*k4/2
    prior = prior/kint

    k1 = qchisq(0.01,kappa1)
    k2 = S1/k1
    k3 = sqrt(k2)
    k4 = 1.2*max(posterior)

    if(plot){
      plot(sigma, prior, type = "l",col = "blue", ylim = c(0,k4), main = expression(paste("Shape of prior and posterior for ", sigma, sep = "")), xlab = expression(sigma),ylab = "Density")
      lines(sigma, posterior, col = "red")
      legend("topleft", lty = 1, lwd = 2, col = c("blue","red"), legend = c("Prior", "Posterior"), bty = "n")
    }
  }else if(kappa<0){
    S0 = 0
    S1 = S0+SST
    kappa1 = kappa+n

    k1 = qchisq(0.001,kappa1)
    k2 = S1/k1
    k3 = sqrt(k2)
    k4 = k3/1000

    sigma = seq(0,k3,length = 1002)[-1]
    sigma.sq = sigma^2

    log.posterior =  -((kappa1-1)/2+1)*log(sigma.sq)-S1/(2*sigma.sq)
    posterior = exp(log.posterior)
    kint = ((2*sum(posterior))-posterior[1]-posterior[1001])*k4/(2*.999)
    posterior = posterior/kint

    log.prior =  -((kappa-1)/2+1)*log(sigma.sq)-S0/(2*sigma.sq)
    prior = exp(log.prior)
    kint = ((2*sum(prior))-prior[1]-prior[1001])*k4/2
    prior = prior/kint

    likelihood = NULL

    k1 = qchisq(0.01,kappa1)
    k2 = S1/k1
    k3 = sqrt(k2)
    k4 = 1.2*max(posterior)

    if(plot){
      plot(sigma, prior, type = "l",col = "blue", xlim = c(0,k3),ylim = c(0,k4),
           main = expression(paste("Shape of prior and posterior for ", sigma, sep = "")),
           xlab = expression(sigma),ylab = "Density")
      lines(sigma, posterior, col = "red")
      legend("topleft", lty = 1, lwd = 2, col = c("blue","red"), legend = c("Prior", "Posterior"), bty = "n")
    }
  }

  cat(paste("S1: ",signif(S1,4)," kappa1 :", signif(kappa1,3),"\n",sep = ""))

  if(cred.int){
    msg = paste0("This argument is deprecated and will not be supported in future releases.",
                 "\nPlease use the quantile function instead.\n")
    warning(msg)
    
    if(kappa1<2)
      cat("Unable to calculate credible interval for sigma if kappa1<= 2\n")
    else{
      sigmahat.post.mean = sqrt(S1/(kappa1-2))
      cat(paste("Estimate of sigma using posterior mean: ",
                signif(sigmahat.post.mean,4),"\n",sep = ""))
    }

    q50 = qchisq(p = 0.5, df = kappa1)
    sigmahat.post.median = sqrt(S1/q50)
    cat(paste("Estimate of sigma using posterior median: ",
              signif(sigmahat.post.median,4),"\n",sep = ""))

    ci = sqrt(S1 / qchisq(p = 1 - c(alpha * 0.5, 1 - alpha * 0.5), df = kappa1))
    ciStr = sprintf("%d%% credible interval for sigma: [%4g, %4g]\n", round(100 * (1 - alpha)), signif(ci[1], 4), signif(ci[2], 4))
    cat(ciStr)
    
    if(plot)
      abline(v = ci, col = "blue", lty = 3)

  }


  results = list(param.x = sigma, prior = prior, likelihood = likelihood,
                 posterior = posterior, 
                 sigma = sigma, # for backwards compat. only
                 S1 = S1, kappa1 = kappa1,
                 mean = ifelse(kappa1 > 2, sqrt(S1 / (kappa1 - 2)), NA),
                 median = sqrt(S1 / qchisq(0.5, kappa1)),
                 var = ifelse(kappa1 > 4, 2 * S1^2 / ((kappa1 - 2)^2 * (kappa1 - 4)), NA),
                 sd = sqrt(ifelse(kappa1 > 4, 2 * S1^2 / ((kappa1 - 2)^2 * (kappa1 - 4)), NA)),
                 cdf = function(y, ...){
                   pchisq(S1 / y^2, df = kappa1, ...)
                 },
                 quantileFun = function(probs, ...){
                   sqrt(S1 / qchisq(p = 1 - probs, df = kappa1, ...))
                 }
                 )
  class(results) = 'Bolstad'
  invisible(results)
}

