#' Bayesian inference on a normal mean with a general continuous prior
#' 
#' Evaluates and plots the posterior density for \eqn{\mu}{mu}, the mean of a
#' normal distribution, with a general continuous prior on \eqn{\mu}{mu}
#' 
#' 
#' @param x a vector of observations from a normal distribution with unknown
#' mean and known std. deviation.
#' @param sigma.x the population std. deviation of the normal distribution
#' @param density distributional form of the prior density can be one of:
#' "normal", "unform", or "user".
#' @param params if density = "normal" then params must contain at least a mean
#' and possible a std. deviation. If a std. deviation is not specified then
#' sigma.x will be used as the std. deviation of the prior. If density =
#' "uniform" then params must contain a minimum and a maximum value for the
#' uniform prior. If a maximum and minimum are not specified then a
#' \eqn{U[0,1]} prior is used
#' @param n.mu the number of possible \eqn{\mu}{mu} values in the prior
#' @param mu a vector of possibilities for the probability of success in a
#' single trial. Must be set if density="user"
#' @param mu.prior the associated prior density. Must be set if density="user"
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced.
#' @return A list will be returned with the following components:
#' \item{likelihood}{the scaled likelihood function for \eqn{\mu}{mu} given
#' \eqn{x} and \eqn{\sigma_x}{sigma.x}} \item{posterior}{the posterior
#' probability of \eqn{\mu}{mu} given \eqn{x} and \eqn{\sigma}{sigma.x}}
#' \item{mu}{the vector of possible \eqn{\mu}{mu} values used in the prior}
#' \item{mu.prior}{the associated probability mass for the values in
#' \eqn{\mu}{mu}}
#' @seealso \code{\link{normdp}} \code{\link{normnp}}
#' @keywords misc
#' @examples
#' 
#' ## generate a sample of 20 observations from a N(-0.5,1) population
#' x = rnorm(20,-0.5,1)
#' 
#' ## find the posterior density with a uniform U[-3,3] prior on mu
#' normgcp(x, 1, params = c(-3, 3))
#' 
#' ## find the posterior density with a non-uniform prior on mu
#' mu = seq(-3, 3, by = 0.1)
#' mu.prior = rep(0, length(mu))
#' mu.prior[mu <= 0] = 1 / 3 + mu[mu <= 0] /9
#' mu.prior[mu > 0] = 1 / 3 - mu[mu > 0] / 9
#' normgcp(x, 1, density = "user", mu = mu, mu.prior = mu.prior)
#' 
#' ## find the CDF for the previous example and plot it
#' ## Note the syntax for sintegral has changed
#' results = normgcp(x,1,density="user",mu=mu,mu.prior=mu.prior)
#' cdf = sintegral(mu,results$posterior,n.pts=length(mu))$cdf
#' plot(cdf,type="l",xlab=expression(mu[0])
#'              ,ylab=expression(Pr(mu<=mu[0])))
#' 
#' ## use the CDF for the previous example to find a 95%
#' ## credible interval for mu. Thanks to John Wilkinson for this simplified code
#' 
#' lcb = cdf$x[with(cdf,which.max(x[y<=0.025]))]
#' ucb = cdf$x[with(cdf,which.max(x[y<=0.975]))]
#' cat(paste("Approximate 95% credible interval : ["
#'            ,round(lcb,4)," ",round(ucb,4),"]\n",sep=""))
#' 
#' ## use the CDF from the previous example to find the posterior mean
#' ## and std. deviation
#' dens = mu*results$posterior
#' post.mean = sintegral(mu,dens)$value
#' 
#' dens = (mu-post.mean)^2*results$posterior
#' post.var = sintegral(mu,dens)$value
#' post.sd = sqrt(post.var)
#' 
#' ## use the mean and std. deviation from the previous example to find
#' ## an approximate 95% credible interval
#' lb = post.mean-qnorm(0.975)*post.sd
#' ub = post.mean+qnorm(0.975)*post.sd
#' 
#' 
#' cat(paste("Approximate 95% credible interval : ["
#'    ,round(lb,4)," ",round(ub,4),"]\n",sep=""))
#' 
#' ## repeat the last example but use the new summary functions for the posterior
#' results = normgcp(x, 1, density = "user", mu = mu, mu.prior = mu.prior)
#' 
#' ## use the cdf function to get the cdf and plot it
#' postCDF = cdf(results) ## note this is a function
#' plot(results$mu, postCDF(results$mu), type="l", xlab = expression(mu[0]),
#'      ylab = expression(Pr(mu <= mu[0])))
#' 
#' ## use the quantile function to get a 95% credible interval
#' ci = quantile(results, c(0.025, 0.975))
#' ci
#' 
#' ## use the mean and sd functions to get the posterior mean and standard deviation
#' postMean = mean(results)
#' postSD = sd(results)
#' postMean
#' postSD
#' 
#' ## use the mean and std. deviation from the previous example to find
#' ## an approximate 95% credible interval
#' ciApprox = postMean + c(-1,1) * qnorm(0.975) * postSD
#' ciApprox
#' 
#' @export normgcp
normgcp = function(x, sigma.x = NULL, density = c("normal", "unform", "user") ,
                   params = NULL, n.mu = 50, mu = NULL,
                   mu.prior = NULL, plot = TRUE){

  ## x - the vector of observations
  ## sigma.x - the population standard deviation
  ## density - distributional form of the prior density
  ## can be one of : normal, unform, or user
  ## by default a continuous uniform prior is used
  ## mu - vector of possible values of the population mean
  ## mu.prior - the associated prior probability mass
  ## ret - if true then the likelihood and posterior are returned as a
  ## list

  mean.x = mean(x)

  if(n.mu < 3)
    stop("Number of prior values of mu must be greater than 2")

  if(is.null(sigma.x)){
    sigma.x = sd(x - mean.x)
    cat(paste("Standard deviation of the residuals :",
              signif(sigma.x,4),"\n", sep = ""))
  }else{
    cat(paste("Known standard deviation :", signif(sigma.x, 4),"\n",sep=""))
  }

  density = match.arg(density)

  if(grepl('^n(orm(al)*)*$', density)){
    if(is.null(params) | length(params) < 1)
      stop("You must supply a mean for a normal prior")
    mx = params[1]

    if(length(params) == 2)  ## user has supplied sd as well
      s.x = params[2]
    else
      s.x = sigma.x

    mu = seq(mx - 3.5 * s.x, mx + 3.5 * s.x, length = n.mu)
    mu.prior = dnorm(mu,mx,s.x)
  }else if(grepl('^u(nif(orm)*)*$', density)){
    if(is.null(params)){
      ## set params to mean+/-3.5sd by default
      params = c(mean.x - 3.5 * sigma.x, mean.x + 3.5 * sigma.x)
    }
    if(length(params)<2)
      stop("You must supply a minimum and a maximum to use a uniform prior")
    minx = params[1]
    maxx = params[2]
    if(maxx <= minx)
      stop("The maximum must be greater than the minimum for a uniform prior")
    mu = seq(minx, maxx, length = n.mu)
    mu.prior = dunif(mu, minx, maxx)
  }else{
    ## user specified prior
    if(is.null(mu) | is.null(mu.prior))
      stop("If you wish to use a non-uniform continuous prior then you must supply a mean vector, mu, and an associated density vector, mu.prior")
    
    if(is.function(mu.prior))
      mu.prior = mu.prior(mu)
  }

  if(any(mu.prior< 0))
    stop("Prior densities must be >=0")

  crude.int = sum(diff(mu) * mu.prior[-1])
  if(round(crude.int, 3) != 1){
    warning("The prior probabilities did not sum to 1, therefore the prior has been normalized")
    mu.prior = mu.prior / crude.int
    print(crude.int)
  }

  n.mu = length(mu)
  mx = mean(x)
  nx = length(x)
  snx = sigma.x^2/nx
  likelihood = exp(-0.5*(mx-mu)^2/snx)

  ## Numerically integrate the denominator
  ## First calculate the height of the function to be integrated

  f.x.mu = likelihood*mu.prior

  ## Now get a linear approximation so that we don't have to worry about
  ## the number of points specified by the user

  ap = approx(mu,f.x.mu,n=513)
  integral = sum(ap$y[2*(1:256)-1]+4*ap$y[2*(1:256)]+ap$y[2*(1:256)+1])
  integral = (ap$x[2]-ap$x[1])*integral/3

  posterior = likelihood*mu.prior/integral

  if(plot){
    plot(mu, posterior, ylim = c(0, 1.1 * max(posterior, mu.prior)), type = "l",
         lty = 1,col="blue",
         xlab = expression(mu), ylab = expression(Probabilty(mu)))
    lines(mu,mu.prior,lty=2,col="red")
  
    legend("topleft", bty = "n", cex = 0.7, 
           lty = 1:2, col = c("blue", "red"),
           legend = c("Posterior", "Prior"))
  }
  
  results = list(name = 'mu', param.x = mu, prior = mu.prior, 
                 likelihood = likelihood, posterior = posterior,
                 mu = mu, mu.prior = mu.prior #for backwards compat. only
  )
  class(results) = 'Bolstad'
  invisible(results)
}
