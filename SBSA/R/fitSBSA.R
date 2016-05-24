#' Conducts sensitivity analysis over a model involving
#' unobserved and poorly measured covariates.
#'
#' The function uses a simplified Bayesian sensitivity analysis
#' algorithm that models the outcome variable \eqn{Y} in terms
#' of exposure \eqn{X} and confounders \eqn{Z=(Z_1,\ldots,Z_p}
#' and \eqn{U=(U_1,\ldots,U_q)}, where \eqn{U}s are unobserved,
#' and \eqn{Z}s are measured imprecisely as \eqn{W}s. (I.e., the
#' observed data is \eqn{(Y, X, W)}.) Parameters of the model
#' are then estimated using MCMC with reparametrizing
#' block-sampling. The estimated parameters are as follows:
#' \itemize{
#' \item \eqn{\tau}: \eqn{(W|Y, U, Z, X) \sim N_p(Z, diag(\tau^2))}
#' \item \eqn{\gamma_x, \gamma_z}: \eqn{(U|X, Z) \sim N(\gamma_x X + \gamma_z' Z)}
#' \item \eqn{\alpha, \beta_u, \beta_z, \sigma}: \eqn{(Y|U, Z, X) \sim N(\alpha_0 + \alpha_x X + \beta_u U + \beta_z' Z, \sigma^2)}
#' }
#'
#' @param y a vector of outcomes
#' @param x a (standardized) vector of exposures
#' @param w a (standardized) matrix of noisy measurements
#' @param a parameter of the prior for magnitude of measurement
#'   error on confounder \eqn{Z_j}
#' @param b parameter of the prior for magnitude of measurement
#'   error on confounder \eqn{Z_j}
#' @param k2 (optional) magnitude of prior uncertainty about
#'   \eqn{(U|X, Z)} regression coefficients
#' @param el2 (optional) residual variance for \eqn{(U|X, Z)} 
#' @param cor.alpha (optional) value of the \eqn{\rho} parameter
#'  of the bivariate normal prior for \eqn{\alpha}
#' @param sd.alpha (optional) value of the \eqn{\sigma} parameter
#'  of the bivariate normal prior for \eqn{\alpha}
#' @param nrep number of MCMC steps
#' @param sampler.jump named vector of standard deviation of
#    jumps for each block sampler with the following named elements:
#'   \itemize{
#'     \item{@code{alpha} jump for block reparametrizing \eqn{\alpha}}
#'     \item{@code{beta.z} jump for block reparametrizing \eqn{\beta_z}}
#'     \item{@code{sigma.sq} (continuous case only) jump for block reparametrizing \eqn{\sigma^2}}
#'     \item{@code{tau.sq} jump for block reparametrizing \eqn{\tau^2}}
#'     \item{@code{beta.u.gamma.x} jump for block reparametrizing \eqn{\beta_u}
#'                          and \eqn{\gamma_z}}
#'     \item{@code{gamma.z} jump for block reparametrizing \eqn{\gamma_z}}
#'   }
#' @param q.steps number of steps in numeric integration of
#' likelihood (only used for binary outcome variables)
#' @param family a character string indicating the assumed
#'   distribution of the outcome. Valid values are \code{"continuous"},
#'   the default, or \code{"binary"}.
#' @return a list with the following elements:
#'   \item{acc}{a vector of counts of how many times each block
#'      sampler successfully made a jump. Vector elements are named by their block,
#'      as in the @code{sampler.jump} argument. }
#'   \item{alpha}{a \eqn{nrep \times\ 2} matrix of the value of \eqn{\alpha} parameter at each MCMC step}
#'   \item{beta.z}{a \eqn{nrep \times\ p} matrix of the value of \eqn{\beta_z} parameter at each MCMC step}
#'   \item{gamma.z}{a \eqn{nrep \times\ p} matrix of the value of \eqn{\gamma_z} parameter at each MCMC step}
#'   \item{tau.sq}{a \eqn{nrep \times\ p} matrix of the value of \eqn{\tau^2} parameter at each MCMC step}
#'   \item{gamma.x}{a vector of the value of \eqn{\gamma_x} parameter at each MCMC step}
#'   \item{beta.u}{a vector of the value of \eqn{\beta_u} parameter at each MCMC step}
#'   \item{sigma.sq}{a vector of the value of \eqn{\sigma^2} parameter at each MCMC step}
#' @title Fitting Simplified Bayesian Sensitivity Models
#' @references
#'   Gustafson, P. and McCandless L. C. and Levy, A. R. and Richardson, S. (2010)
#'   \emph{Simplified Bayesian Sensitivity Analysis for Mismeasured
#'         and Unobserved Confounders.}
#'   Biometrics, 66(4):1129--1137.
#'   DOI: 10.1111/j.1541-0420.2009.01377.x
#' @keywords TODO
#' @examples
#' ### simulated data example
#' n <- 1000
#' 
#' ### exposure and true confounders equi-correlated with corr=.6
#' tmp <- sqrt(.6)*matrix(rnorm(n),n,5) +
#'        sqrt(1-.6)*matrix(rnorm(n*5),n,5)
#' x <- tmp[,1]
#' z <- tmp[,2:5]
#' 
#' ### true outcome relationship
#' y <- rnorm(n, x + z%*%rep(.5,4), .5)
#' 
#' 
#' ### first two confounders are poorly measured, ICC=.7, .85
#' ### third is correctly measured, fourth is unobserved
#' w <- z[,1:3]
#' w[,1] <- w[,1] + rnorm(n, sd=sqrt(1/.7-1))
#' w[,2] <- w[,2] + rnorm(n, sd=sqrt(1/.85-1))
#' 
#' ### fitSBSA expects standardized exposure, noisy confounders
#' x.sdz <- (x-mean(x))/sqrt(var(x))
#' w.sdz <- apply(w, 2, function(x) {(x-mean(x)) / sqrt(var(x))} )
#' 
#' ### prior information: ICC very likely above .6, mode at .8
#' ### via Beta(5,21) distribution
#' fit <- fitSBSA(y, x.sdz, w.sdz, a=5, b=21, nrep=20000,
#'                sampler.jump=c(alpha=.02, beta.z=.03,
#'                               sigma.sq=.05, tau.sq=.004,
#'                               beta.u.gamma.x=.4, gamma.z=.5))
#' 
#' ### check MCMC behaviour
#' print(fit$acc)
#' plot(fit$alpha[,2], pch=20)
#' 
#' ### inference on target parameter in original scale
#' trgt <- fit$alpha[10001:20000,2]/sqrt(var(x))
#' print(c(mean(trgt), sqrt(var(trgt))))
fitSBSA <- function(y, x, w, a, b,
                   k2=NULL, el2=NULL, cor.alpha=0, sd.alpha=1e6, 
                   nrep=5000, sampler.jump=c(alpha=.15, beta.z=.1,
                                sigma.sq=.5, tau.sq=.05,
                                beta.u.gamma.x=.3, gamma.z=.15),
                   q.steps=25,
                   family=c('continuous', 'binary')) {
  family <- match.arg(family)

  # Check that all the expected sampler.jumps are present
  required.samplers <- switch(family,
                              continuous=c('alpha', 'beta.z',
                                'sigma.sq', 'tau.sq',
                                'beta.u.gamma.x', 'gamma.z'),
                              binary=c('alpha', 'beta.z', 'tau.sq',
                                'beta.u.gamma.x', 'gamma.z'))
  if (length(missing.samplers <- setdiff(required.samplers,
                                         names(sampler.jump)))) {
    stop('sampler.jump argument is missing required elements: ',
         paste(missing.samplers, ''))
  }
  
  if (length(sd.alpha) == 1) {
    sd.alpha.0 <- sd.alpha.x <- sd.alpha
  }
  else if (length(sd.alpha) == 2) {
    sd.alpha.0 <- sd.alpha[1]
    sd.alpha.x <- sd.alpha[1]
  }
  else {
    stop("'sd.alpha' must be a one or two-element vector")
  }

  p <- dim(w)[2]

  # Check that W is standardized by columns
  if (!isTRUE(all.equal(apply(w, 2, mean), rep(0, p))) ||
      !isTRUE(all.equal(apply(w, 2, var), rep(1, p)))) {
    warning('The W matrix should be standardized!!')
  }
  
  if (length(a) == 1) {
    a <- rep(a, p)
  }
  else if (length(a) != p) {
    stop("'a' must be a scalar or a vector of same length as 'y'")
  }
  
  if (length(b) == 1) {
    b <- rep(b, p)
  }
  else if (length(b) != p) {
    stop("'b' must be a scalar or a vector of same length as 'y'")
  }
  
  Sigma <- var(cbind(x,w))  
  M <- Sigma[-1,-1]-(1/Sigma[1,1])*Sigma[1,-1]%*%t(Sigma[1,-1])
  mu <- Sigma[1,-1]/Sigma[1,1]
  
  if ( is.null(k2) && is.null(el2) ) {
    k2 <- mean((a+b)/b) * (sum(Sigma^2)-sum(diag(Sigma^2))) / ((p+1)^2-(p+1))
    el2 <- max(mean(b/(a+b)) - (p+1)*k2, 1)
  }
  if (el2 <= 0) {
    stop(paste("'el2' argument must be positive:", el2))
  }
  if (k2 <= 0) {
    stop(paste("'k2' argument must be positive:", k2))
  }

  sampler.jump <- as.list(sampler.jump[required.samplers])
  for (block in c('beta.z', 'tau.sq', 'gamma.z')) {
    sampler.jump[[block]] <- rep(sampler.jump[[block]], length.out=p)
  }
  sampler.jump[['alpha']] <- rep(sampler.jump[['alpha']], 2)
  
  ## TODO: name sampler jump elements
  switch(family,
         continuous=.Call('fitbsa', y, x, w,
                          a, b, k2, el2, nrep,
                          sampler.jump,
                          Sigma, M, mu, cor.alpha,
                          sd.alpha.0, sd.alpha.x,
                          PACKAGE = 'SBSA'),
         
         binary=.Call('fitbsa_binary', y, x, w,
                      a, b, k2, el2, nrep,
                      sampler.jump,
                      Sigma, M, mu,
                      q.steps, qnorm((1:q.steps)/(q.steps+1)),
                      cor.alpha,
                      sd.alpha.0, sd.alpha.x,
                      PACKAGE = 'SBSA'))
}
