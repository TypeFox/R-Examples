#' Poisson sampling with a general continuous prior
#' 
#' Evaluates and plots the posterior density for \eqn{\mu}{mu}, the mean rate
#' of occurance of an event or objects, with Poisson sampling and a general
#' continuous prior on \eqn{\mu}{mu}
#' 
#' 
#' @param y A random sample of one or more observations from a Poisson
#' distribution
#' @param density may be one of "gamma", "normal", or "user"
#' @param params if density is one of the parameteric forms then then a vector
#' of parameters must be supplied.  gamma: a0,b0 normal: mean,sd
#' @param n.mu the number of possible \eqn{\mu}{mu} values in the prior. This
#' number must be greater than or equal to 100. It is ignored when
#' density="user".
#' @param mu a vector of possibilities for the mean of a Poisson distribution.
#' This must be set if density="user".
#' @param mu.prior the associated prior density. This must be set if
#' density="user".
#' @param print.sum.stat if set to TRUE then the posterior mean, posterior
#' variance, and a credible interval for the mean are printed. The width of the
#' credible interval is controlled by the parameter alpha.
#' @param alpha The width of the credible interval is controlled by the
#' parameter alpha.
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced.
#' @param suppressOutput if \code{TRUE} then none of the output is printed to
#' console
#' @return A list will be returned with the following components: \item{mu}{the
#' vector of possible \eqn{\mu}{mu} values used in the prior}
#' \item{mu.prior}{the associated probability mass for the values in
#' \eqn{\mu}{mu}} \item{likelihood}{the scaled likelihood function for
#' \eqn{\mu}{mu} given \eqn{y}} \item{posterior}{the posterior probability of
#' \eqn{\mu}{mu} given \eqn{y}}
#' @seealso \code{\link{poisdp}} \code{\link{poisgamp}}
#' @keywords misc
#' @examples
#' 
#' ## Our data is random sample is 3, 4, 3, 0, 1. We will try a normal
#' ## prior with a mean of 2 and a standard deviation of 0.5.
#' y = c(3,4,3,0,1)
#' poisgcp(y,density="normal",params=c(2,0.5))
#' 
#' ## The same data as above, but with a gamma(6,8) prior
#' y = c(3,4,3,0,1)
#' poisgcp(y,density="gamma",params=c(6,8))
#' 
#' ## The same data as above, but a user specified continuous prior.
#' ## We will use print.sum.stat to get a 99% credible interval for mu.
#' y = c(3,4,3,0,1)
#' mu = seq(0,8,by=0.001)
#' mu.prior = c(seq(0,2,by=0.001),rep(2,1999),seq(2,0,by=-0.0005))/10
#' poisgcp(y,"user",mu=mu,mu.prior=mu.prior,print.sum.stat=TRUE,alpha=0.01)
#' 
#' ## find the posterior CDF using the results from the previous example
#' ## and Simpson's rule. Note that the syntax of sintegral has changed.
#' results = poisgcp(y,"user",mu=mu,mu.prior=mu.prior)
#' cdf = sintegral(mu,results$posterior,n.pts=length(mu))$cdf
#' plot(cdf,type="l",xlab=expression(mu[0])
#' 	,ylab=expression(Pr(mu<=mu[0])))
#' 
#' ## use the cdf to find the 95% credible region.
#' lcb = cdf$x[with(cdf,which.max(x[y<=0.025]))]
#' ucb = cdf$x[with(cdf,which.max(x[y<=0.975]))]
#' cat(paste("Approximate 95% credible interval : ["
#' 	,round(lcb,4)," ",round(ucb,4),"]\n",sep=""))
#' 
#' ## find the posterior mean, variance and std. deviation
#' ## using Simpson's rule and the output from the previous example
#' dens = mu*results$posterior # calculate mu*f(mu | x, n)
#' post.mean = sintegral(mu,dens)$value
#' 
#' dens = (mu-post.mean)^2*results$posterior
#' post.var = sintegral(mu,dens)$value
#' post.sd = sqrt(post.var)
#' 
#' # calculate an approximate 95% credible region using the posterior mean and
#' # std. deviation
#' lb = post.mean-qnorm(0.975)*post.sd
#' ub = post.mean+qnorm(0.975)*post.sd
#' 
#' cat(paste("Approximate 95% credible interval : ["
#' 	,round(lb,4)," ",round(ub,4),"]\n",sep=""))
#' 
#' @export poisgcp
poisgcp = function(y, density = c("normal", "gamma", "user"), params = c(0, 1), n.mu = 100,
                  mu = NULL, mu.prior = NULL,
                  print.sum.stat = FALSE, alpha = 0.05, plot = TRUE,
                  suppressOutput = FALSE){
  n = length(y)
  y.sum = sum(y)
  y.bar = mean(y)

  if(is.null(y))
    stop("Error: y has no data")

  if(any(y < 0))
    stop("Error: data contains negative values")

  if(n.mu<100)
    stop("Error: there must be at least 100 points in the prior")

  density = match.arg(density)
  
  if(density=="user"){
    if(is.null(mu) || is.null(mu.prior))
      stop("Error: a vector of possibilities (mu) and associated densities must be specified for a user prior")

    if(length(mu)!=length(mu.prior))
      stop("Error: There must be an equal number of values in mu and mu prior")

  }else if(density=="normal"){
    if(length(params)!=2)
      stop("Error: A mean and a std. deviation must be specified for a normal prior")
    mx = params[1]
    sx = params[2]

    if(sx<=0)
      stop("Error: the std. deviation of a normal prior must be greater than zero")

    lb = mx-3.5*sx
    ub = mx+3.5*sx

    if(lb<0)
      { cat("The normal prior has negative values.")
        cat("Whilst this is true for all normal distributions, you can sneak\n")
        cat("around it by using a large positive mean relative to the \n")
        cat("std. deviation.\n")
        stop("Error")
      }

    mu = seq(lb,ub,length=n.mu)
    mu.prior = dnorm(mu,mx,sx)
  }else if(density=="gamma"){
    if(length(params)!=2)
      stop("Error: there must be two parameters, a0 and b0 for a gamma prior")
    if(sum(params<0)>0)
      stop("Error: the parameters of a gamma prior must be positive")

    a0 = params[1]
    b0 = params[2]
    gamma.bds = qgamma(c(0.005,0.995),a0,b0)
    mu = seq(gamma.bds[1],gamma.bds[2],length=n.mu)
    mu.prior = dgamma(mu,a0,b0)
  }else{
    stop(paste("Error: unrecognized density: ",density,". The options are normal, gamma or user."))
  }

  if(sum(mu<0)>0)
    stop("Error: mu cannot contain negative values")

  if(!suppressOutput){
    cat("Summary statistics for data\n")
    cat("---------------------------\n")
    cat(paste("Number of observations:\t", n,"\n"))
    cat(paste("Sum of observations:\t", y.sum,"\n"))
  }

  log.lik = y.sum*log(mu)-n*mu
  likelihood = exp(log.lik)
  fx.joint = approxfun(mu,mu.prior*likelihood)
  normalizing.constant = integrate(fx.joint,min(mu),max(mu))$value
  posterior = likelihood*mu.prior/normalizing.constant

  if(print.sum.stat){
    fx.posterior = approxfun(mu,posterior)
    x.fx = approxfun(mu,posterior*mu)
    posterior.mean = integrate(x.fx,min(mu),max(mu))$value
    xmusq.fx = approxfun(mu,(mu-posterior.mean)^2*posterior)
    posterior.var = integrate(xmusq.fx,min(mu),max(mu))$value
    cat("\nPosterior distribution summary statistics\n")
    cat("-----------------------------------------\n")
    cat(paste("Post. mean:\t", round(posterior.mean,3), "\n"))
    cat(paste("Post. var.:\t", round(posterior.var,4), "\n"))

    mu.int = seq(min(mu),max(mu),length=256)
    f.mu = fx.posterior(mu.int)
    suppressMessages({cdf = sintegral(mu.int,f.mu)$cdf;
                      fx.posterior.invcdf = approxfun(cdf$y,cdf$x)})
    lb = fx.posterior.invcdf(alpha/2)
    ub = fx.posterior.invcdf(1-alpha/2)
    cat(paste(round(100*(1-alpha)),"% cred. int.: ["
              , round(lb,3), ",", round(ub,3),"]\n\n"))
  }

  if(plot){
    y.max = max(mu.prior, posterior)
    plot(mu,mu.prior,ylim=c(0,1.1*y.max),xlab=expression(mu)
         ,ylab="Density",
         ,main="Shape of continuous prior and posterior for Poisson mean"
         ,type="l",lty=2,col="red")
    lines(mu,posterior,lty=3,col="blue")
    legend("topleft", bty = "n", lty = 2:3, col=c("red", "blue"),
           legend=c("Prior","Posterior"), cex = 0.7)
  }
  
  results = list(name =  'mu', param.x = mu, prior = mu.prior, likelihood = likelihood, posterior = posterior,
                 mu = mu, mu.prior = mu.prior # for backwards compatibility only
                 )
  class(results) = 'Bolstad'
  invisible(results)
}



