#' Poisson sampling with a discrete prior
#' 
#' Evaluates and plots the posterior density for \eqn{\mu}{mu}, the mean rate
#' of occurance in a Poisson process and a discrete prior on \eqn{\mu}{mu}
#' 
#' 
#' @param y.obs a random sample from a Poisson distribution.
#' @param mu a vector of possibilities for the mean rate of occurance of an
#' event over a finite period of space or time.
#' @param mu.prior the associated prior probability mass.
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced.
#' @return A list will be returned with the following components:
#' 
#' \item{likelihood}{the scaled likelihood function for \eqn{\mu}{mu} given
#' \eqn{y_{obs}}{y.obs}} \item{posterior}{the posterior probability of
#' \eqn{\mu}{mu} given \eqn{y_{obs}}{y.obs}} \item{mu}{the vector of possible
#' \eqn{\mu}{mu} values used in the prior} \item{mu.prior}{the associated
#' probability mass for the values in \eqn{\mu}{mu}}
#' @seealso \code{\link{poisgamp}} \code{\link{poisgcp}}
#' @keywords misc
#' @examples
#' 
#' ## simplest call with an observation of 4 and a uniform prior on the
#' ## values mu = 1,2,3
#' poisdp(4,1:3,c(1,1,1)/3)
#' 
#' ##  Same as the previous example but a non-uniform discrete prior
#' mu = 1:3
#' mu.prior = c(0.3,0.4,0.3)
#' poisdp(4,mu=mu,mu.prior=mu.prior)
#' 
#' ##  Same as the previous example but a non-uniform discrete prior
#' mu = seq(0.5,9.5,by=0.05)
#' mu.prior = runif(length(mu))
#' mu.prior = sort(mu.prior/sum(mu.prior))
#' poisdp(4,mu=mu,mu.prior=mu.prior)
#' 
#' ## A random sample of 50 observations from a Poisson distribution with
#' ## parameter mu = 3 and  non-uniform prior
#' y.obs = rpois(50,3)
#' mu = c(1:5)
#' mu.prior = c(0.1,0.1,0.05,0.25,0.5)
#' results = poisdp(y.obs, mu, mu.prior)
#' 
#' ##  Same as the previous example but a non-uniform discrete prior
#' mu = seq(0.5,5.5,by=0.05)
#' mu.prior = runif(length(mu))
#' mu.prior = sort(mu.prior/sum(mu.prior))
#' y.obs = rpois(50,3)
#' poisdp(y.obs,mu=mu,mu.prior=mu.prior)
#' 
#' 
#' @export poisdp
poisdp = function(y.obs, mu, mu.prior, plot = TRUE){
  if(length(y.obs) == 0 || is.null(y.obs))
    stop("Error: y.obs must contain at least one value")

  if(any(y.obs < 0))
    stop("Error: y.obs cannot contain negative values")

  if(length(mu) != length(mu.prior))
    stop("Error: the lengths of mu and mu.prior are unequal.\nThere must be a corresponding probability for each value of mu")

  if(sum(mu<=0)>0)
    stop("Error: the values of the rate paramter mu, must be greater than zero")

  if(sum(mu.prior<0)>0 || sum(mu.prior>1)>0)
    stop("Error: prior probabilities must be between zero and one")

  if(sum(mu.prior)!=1){
    warning("Warning: the prior does not sum to 1. The prior has been rescaled")
    mu.prior = mu.prior/sum(mu.prior)
  }

  n = length(y.obs)
  if(n == 1){
    k = y.obs

    m = length(mu)
    cat("Prior\n")
    cat("-----\n")
    prior.matrix = cbind(mu,mu.prior)
    colnames(prior.matrix) = c("mu","Pr(mu)")
    print(prior.matrix)
    k1 = 0.9995
    k2 = mu[m]
    cat(paste("\nk1:\t",k1,"\nk2:\t",k2,"\n\n"))

    n = qpois(k1,k2)
    y1 = 0:n

    k1 = mu
    k1[k1==0] = 1e-9

    f.cond = matrix(0,nrow=m,ncol=n+1)
    for(i in 1:m)
      f.cond[i,] = dpois(y1,k1[i])

    rownames(f.cond) = mu
    colnames(f.cond) = y1

    cat("Conditional probability of y1 given mu\n")
    cat("-------------------------------------\n")
    print(f.cond)
    cat("\n\n")

    matrix.prior = diag(mu.prior)
    f.joint = matrix.prior%*%f.cond

    cat("Joint probability of y1 and mu\n")
    cat("-------------------------------\n")
    print(f.joint)
    cat("\n\n")

    f.marg = apply(f.joint,2,sum)

    cat("Marginal probability of y1\n")
    cat("-------------------------\n")
    print(f.marg)
    cat("\n\n")

    ## extract the column of f.cond corresponding to y.obs
    likelihood = f.cond[,y.obs+1]
    posterior = likelihood*mu.prior
    posterior = posterior/sum(posterior)

    results = cbind(mu,mu.prior,likelihood,posterior)
    colnames(results) = c("Mu","Prior","Likelihood","Posterior")
    print(results)
  }else{
    m = length(mu)

    likelihood = rep(0,m)
    for(i in 1:m)
      likelihood[i] = exp(sum(dpois(y.obs,max(mu[i],1e-9),log=TRUE)))

    posterior = likelihood*mu.prior
    posterior = posterior/sum(posterior)

    results = cbind(mu,mu.prior,likelihood,posterior)
    colnames(results) = c("Mu","Prior","Likelihood","Posterior")
    print(results)
  }

  if(plot){
    plot.data = rbind(mu.prior,posterior)
    if(length(mu.prior)<=10){
      colnames(plot.data) = mu
      y.max = max(mu.prior,posterior)
      midpoints = barplot(plot.data,beside=TRUE,col=c("red","blue")
                         ,xlab=expression(mu)
                         ,ylab=expression(paste("Pr(",mu,"|","y)"))
                         ,ylim=c(0,y.max*1.1)
                         ,main=expression(
                             paste("Prior and posterior probability for ", mu
                                   ," given the data y")))
      legend("topleft", cex = 0.7, bty = "n", 
             legend=c("Prior","Posterior"),
             fill=c("red","blue"))
      box()
    }else{
      y.max = max(mu.prior,posterior)
      plot(mu,mu.prior,type="l",lty=2,col="red"
           ,xlab=expression(mu)
           ,ylab=expression(paste("Pr(",mu,"|","y)"))
           ,ylim=c(0,y.max*1.1)
           ,main=expression(paste("Prior and posterior probability for ", mu
               ," given the data y")))
      lines(mu,posterior,lty=1,col="blue")
      legend("topleft", cex = 0.7, bty = "n", 
             lty = c(2,1), col = c("red"," blue"),
             legend=c("Prior","Posterior"))
  
    }
  }
  mx = sum(mu * posterior)
  vx = sum((mu - mx)^2 * posterior)
  
  results= list(name = 'mu', param.x = mu, prior = mu.prior, likelihood = likelihood, posterior = posterior,
                mean = mx, var = vx, 
                cdf = function(m){cumDistFun(m, mu, posterior)},
                quantileFun = function(probs, ...){qFun(probs, mu, posterior)},
                mu = mu, mu.prior = mu.prior #for backwards compat. only
                )
  class(results) = 'Bolstad'
  invisible(results)
}



