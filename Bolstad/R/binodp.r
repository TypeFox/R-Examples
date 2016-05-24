#' Binomial sampling with a discrete prior
#' 
#' Evaluates and plots the posterior density for \eqn{\pi}{pi}, the probability
#' of a success in a Bernoulli trial, with binomial sampling and a discrete
#' prior on \eqn{\pi}{pi}
#' 
#' 
#' @param x the number of observed successes in the binomial experiment.
#' @param n the number of trials in the binomial experiment.
#' @param pi a vector of possibilities for the probability of success in a
#' single trial. if \code{pi} is \code{NULL} then a discrete uniform prior for
#' \eqn{\pi}{pi} will be used.
#' @param pi.prior the associated prior probability mass.
#' @param n.pi the number of possible \eqn{\pi}{pi} values in the prior
#' @param plot if \code{TRUE} then a plot showing the prior and the posterior
#' will be produced
#' @param suppressOutput if \code{TRUE} then none of the output is printed to
#' console
#' @return A list will be returned with the following components: \item{pi}{the
#' vector of possible \eqn{\pi}{pi} values used in the prior}
#' \item{pi.prior}{the associated probability mass for the values in
#' \eqn{\pi}{pi}} \item{likelihood}{the scaled likelihood function for
#' \eqn{\pi}{pi} given \eqn{x} and \eqn{n}} \item{posterior}{the posterior
#' probability of \eqn{\pi}{pi} given \eqn{x} and \eqn{n}} \item{f.cond}{the
#' conditional distribution of \eqn{x} given \eqn{\pi}{pi} and \eqn{n}}
#' \item{f.joint}{the joint distribution of \eqn{x} and \eqn{\pi}{pi} given
#' \eqn{n}} \item{f.marg}{the marginal distribution of \eqn{x}}
#' @seealso \code{\link{binobp}} \code{\link{binogcp}}
#' @keywords misc
#' @examples
#' 
#' ## simplest call with 6 successes observed in 8 trials and a uniform prior
#' binodp(6,8)
#' 
#' ## same as previous example but with more possibilities for pi
#' binodp(6, 8, n.pi = 100)
#' 
#' ## 6 successes, 8 trials and a non-uniform discrete prior
#' pi = seq(0, 1, by = 0.01)
#' pi.prior = runif(101)
#' pi.prior = sort(pi.prior / sum(pi.prior))
#' binodp(6, 8, pi, pi.prior)
#' 
#' ## 5 successes, 6 trials, non-uniform prior
#' pi = c(0.3, 0.4, 0.5)
#' pi.prior = c(0.2, 0.3, 0.5)
#' results = binodp(5, 6, pi, pi.prior)
#' 
#' ## plot the results from the previous example using a side-by-side barplot
#' results.matrix = rbind(results$pi.prior,results$posterior)
#' colnames(results.matrix) = pi
#' barplot(results.matrix, col = c("red", "blue"), beside = TRUE,
#' 	      xlab = expression(pi), ylab=expression(Probability(pi)))
#' box()
#' legend("topleft", bty = "n", cex = 0.7, 
#'        legend = c("Prior", "Posterior"), fill = c("red", "blue"))
#' 
#' @export binodp
binodp = function(x, n, pi = NULL, pi.prior = NULL, n.pi = 10, 
                  plot = TRUE, suppressOutput = FALSE){

  ## n - the number of trials in the binomial
  ## x - the number of observed successes
  ## pi - the probability of success
  ## pi.prior - the associated prior probability mass
  ## ret - if true then the likelihood and posterior are returned as a
  ## list

  if(x > n)
    stop("The number of observed successes (x) must be smaller than the number of trials (n)")
  
  if(n.pi < 3)
    stop("Number of prior values of pi must be greater than 2")

  if(is.null(pi) | is.null(pi.prior)){
    pi = seq(0,1, length = n.pi)
    pi.prior = rep(1 / n.pi, n.pi)
  }

  if(any(pi < 0) | any(pi > 1))   ## check that probabilities lie on [0,1]
    stop("Values of pi must be between 0 and 1 inclusive")

  if(any(pi.prior < 0) | any(pi.prior > 1))
    stop("Prior probabilities must be between 0 and 1 inclusive")

  if(round(sum(pi.prior), 7) != 1){
    warning("The prior probabilities did not sum to 1, therefore the prior has been normalized")
    pi.prior = pi.prior / sum(pi.prior)
  }
  
  ## make sure possible values are in ascending order
  o = order(pi)
  pi = pi[o]
  pi.prior = pi.prior[o]

  n.pi = length(pi)

  likelihood = dbinom(x,n,pi)
  lp = likelihood * pi.prior
  posterior = lp / sum(lp)

  if(plot){
    plot(pi, posterior,ylim = c(0, 1.1 * max(posterior, pi.prior)),pch=20
       ,col="blue",
       xlab = expression(pi), ylab = expression(Probabilty(pi)))
    points(pi,pi.prior,pch=20,col="red")

    legend("topleft", bty = "n", fill = c("blue", "red"),
           legend = c("Posterior", "Prior"), cex = 0.7)
  }
  ## calculate the Conditional distribution

  f.cond = matrix(0,nrow=n.pi,ncol=n+1)
  rownames(f.cond) = as.character(round(pi,3))
  colnames(f.cond) = as.character(0:n)

  for(i in 1:n.pi)
    f.cond[i,] = dbinom(0:n,n,pi[i])

 
  ## caculate the joint distribution of pi and x given n

  f.joint = diag(pi.prior)%*%f.cond
 
  ## calculate the marginal distribtion

  f.marg = matrix(1,nrow=1,ncol=n.pi)%*%f.joint
  
  
  if(!suppressOutput){
    cat("Conditional distribution of x given pi and  n:\n\n")
    print(round(f.cond,4))
  
    cat("\nJoint distribution:\n\n")
    print(round(f.joint,4))
    
  
    cat("\nMarginal distribution of x:\n\n")
    print(round(f.marg,4))
    cat("\n\n")
  
    ## finally display the prior, likelihood, and posterior
  
    results = cbind(pi.prior,likelihood,posterior)
    rownames(results) = as.character(round(pi,3))
    colnames(results) = c("Prior","Likelihood","Posterior")
  
    print(results)
  }
  
  mx = sum(pi * posterior)
  vx = sum((pi  - mx)^2 * posterior)

  results = list(name = 'pi', param.x = pi, prior = pi.prior, likelihood = likelihood,
                posterior = posterior,
                mean = mx,
                var = vx, 
                cdf = function(X, ...){cumDistFun(X, pi, posterior)},
                quantileFun = function(probs, ...){qFun(probs, pi, posterior)},
                pi = pi, pi.prior = pi.prior, ## this duplication is for backward compatibility
                f.cond = f.cond, f.joint = f.joint, f.marg = f.marg)
  class(results) = 'Bolstad'
  invisible(results)
}
