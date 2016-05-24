#' Binomial sampling with a general continuous prior
#' 
#' Evaluates and plots the posterior density for \eqn{\pi}{pi}, the probability
#' of a success in a Bernoulli trial, with binomial sampling and a general
#' continuous prior on \eqn{\pi}{pi}
#' 
#' 
#' @param x the number of observed successes in the binomial experiment.
#' @param n the number of trials in the binomial experiment.
#' @param density may be one of "beta", "exp", "normal", "student", "uniform"
#' or "user"
#' @param params if density is one of the parameteric forms then then a vector
#' of parameters must be supplied.  beta: a, b exp: rate normal: mean, sd
#' uniform: min, max
#' @param n.pi the number of possible \eqn{\pi}{pi} values in the prior
#' @param pi a vector of possibilities for the probability of success in a
#' single trial. This must be set if density = "user".
#' @param pi.prior the associated prior probability mass. This must be set if
#' density = "user".
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced.
#' @return A list will be returned with the following components:
#' \item{likelihood}{the scaled likelihood function for \eqn{\pi}{pi} given
#' \eqn{x} and \eqn{n}} \item{posterior}{the posterior probability of
#' \eqn{\pi}{pi} given \eqn{x} and \eqn{n}} \item{pi}{the vector of possible
#' \eqn{\pi}{pi} values used in the prior} \item{pi.prior}{the associated
#' probability mass for the values in \eqn{\pi}{pi}}
#' @seealso \code{\link{binobp}} \code{\link{binodp}}
#' @keywords misc
#' @examples
#' 
#' ## simplest call with 6 successes observed in 8 trials and a continuous
#' ## uniform prior
#' binogcp(6, 8)
#' 
#' ## 6 successes, 8 trials and a Beta(2, 2) prior
#' binogcp(6, 8,density = "beta", params = c(2, 2))
#' 
#' ## 5 successes, 10 trials and a N(0.5, 0.25) prior
#' binogcp(5, 10, density = "normal", params = c(0.5, 0.25))
#' 
#' ## 4 successes, 12 trials with a user specified triangular continuous prior
#' pi = seq(0, 1,by = 0.001)
#' pi.prior = rep(0, length(pi))
#' priorFun = createPrior(x = c(0, 0.5, 1), wt = c(0, 2, 0))
#' pi.prior = priorFun(pi)
#' results = binogcp(4, 12, "user", pi = pi, pi.prior = pi.prior)
#' 
#' ## find the posterior CDF using the previous example and Simpson's rule
#' myCdf = cdf(results)
#' plot(myCdf, type = "l", xlab = expression(pi[0]),
#' 	   ylab = expression(Pr(pi <= pi[0])))
#' 
#' ## use the quantile function to find the 95% credible region.
#' qtls = quantile(results, probs = c(0.025, 0.975))
#' cat(paste("Approximate 95% credible interval : ["
#' 	, round(qtls[1], 4), " ", round(qtls, 4), "]\n", sep = ""))
#' 
#' ## find the posterior mean, variance and std. deviation
#' ## using the output from the previous example
#' post.mean = mean(results)
#' post.var = var(results)
#' post.sd = sd(results)
#' 
#' # calculate an approximate 95% credible region using the posterior mean and
#' # std. deviation
#' lb = post.mean - qnorm(0.975) * post.sd
#' ub = post.mean + qnorm(0.975) * post.sd
#' 
#' cat(paste("Approximate 95% credible interval : ["
#' 	, round(lb, 4), " ", round(ub, 4), "]\n", sep = ""))
#' 
#' @export binogcp
binogcp = function(x, n, density = c("uniform", "beta", "exp", "normal",  "user"), 
                   params = c(0, 1), n.pi = 1000, 
                   pi = NULL, pi.prior = NULL, plot = TRUE){

  ## n - the number of trials in the binomial
  ## x - the number of observed successes
  ## density - may be one of "exp", "normal", "uniform" or "user"
  ## params - if the density is not "user" then a vector of parameters
  ## must be supplied.
  ##	exp:		rate
  ##	normal: 	mean, sd
  ##	uniform: 	min, max
  ## n.pi - the number of points to divide the [0, 1] interval into

  ## pi and pi.prior are only specified if density == "user"
  ## pi - the probability of success
  ## pi.prior - the associated prior probability mass
  ## plot - if true then the likelihood and posterior are returned as a
  ## list

  if(x > n)
    stop("The number of observed successes (x) must be smaller than the number of trials")
  if(n.pi < 100)
    stop("Number of prior values of pi must be greater than 100")

  if(is.null(pi) || is.null(pi.prior))
    pi = ppoints(n.pi)
  else{
    if(length(pi) != length(pi.prior))
      stop("pi and pi.prior must have same length")

    if(any(pi < 0))    ## check that the density values are greater than 0
      stop("Values of pi must be >= 0")
  }
  
  density = match.arg(density)

  if(density == "beta"){
    if(length(params) < 2){
      warning("Beta prior requires two shape parameters. Default value Beta(1, 1) = Uniform is being used")
      a = 1
      b = 1
    }else{
      if(params[1] <= 0 | params[2] <= 0)
        stop("Beta prior shape parameters must be greater than zero")
      a = params[1]
      b = params[2]
    }
    pi.prior = dbeta(pi, a,b)
  }else	if(density == "exp"){
    if(params[1] <= 0){
      stop("Parameter for exponential density must be greater than zero")
    }else{
      rate = params[1]
      pi.prior = dexp(pi, rate)
    }
  }else if(density == "normal"){
    if(length(params) < 2)
      stop("Normal prior requires a mean and std. deviation")
    else{
      mx = params[1]
      sx = params[2]
      if(sx <= 0)
        stop("Std. deviation for normal prior must be greater than zero")
      pi.prior = dnorm(pi, mx, sx)
    }
  }else if(density == "uniform"){
    if(length(params) < 2)
      stop("Uniform prior requires a minimum and a maximum")
    else{
      minx = params[1]
      maxx = params[2]

      if(maxx <= minx)
        stop("Maximum must be greater than minimum for a uniform prior")
      pi.prior = dunif(pi, minx, maxx)
    }
  }else if (density != "user"){
    stop(paste("Unrecognized density :", density))
  }

  likelihood = (pi^x) * ((1 - pi)^(n - x))

  ## Numerically integrate the denominator
  ## First calculate the height of the function to be integrated

  f.x.pi = likelihood * pi.prior

  ## Now get a linear approximation so that we don't have to worry about
  ## the number of points specified by the user

  ap = approx(pi, f.x.pi, n = 513)
  integral = sum(ap$y[2 * (1:256) - 1] + 4 * ap$y[2 * (1:256)] + ap$y[2 * (1:256) + 1])
  integral = (ap$x[2] - ap$x[1]) * integral / 3

  posterior = likelihood * pi.prior / integral

  if(plot){
    plot(pi, posterior, ylim = c(0, 1.1 * max(posterior, pi.prior)), lty = 1, type = "l", col = "blue", 
         xlab = expression(pi), ylab = "Density")
    lines(pi, pi.prior, lty = 2, col = "red")
  
    left = min(pi) + diff(range(pi)) * 0.05
    legend("topleft", bty = "n", lty = 1:2, col = c("blue", "red"), 
           legend = c("Posterior", "Prior"), cex = 0.7)
  }
  
  results = list(name = 'pi', param.x = pi, prior = pi.prior, likelihood = likelihood, posterior = posterior, 
                 pi = pi, pi.prior = pi.prior)
  class(results) = 'Bolstad'
  invisible(results)
}
