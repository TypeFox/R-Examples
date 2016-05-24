# From SamplerCompare, (c) 2010 Madeleine Thompson

# distributions.R contains example distributions and functions for
# generating example distributions.

#### GAUSSIANS ####

# See ?make.gaussian for more information.

make.gaussian <- function(mean, sigma=NULL, rho=NULL) {
  ndim <- length(mean)

  # If sigma and rho are unset, set sigma to the identity matrix.
  # If sigma isn't but rho is, make sure rho would lead to a positive
  # definite matrix, then use it to make sigma.

  if (is.null(sigma) && is.null(rho)) {
    sigma <- diag(1, nrow=ndim)
  } else if (is.null(sigma)) {
    stopifnot(rho<1 && rho>-1/(ndim-1))
    sigma <- array(rho, c(ndim, ndim))
    diag(sigma) <- 1
  }

  stopifnot(ndim==dim(sigma)[1] &&
            ndim==dim(sigma)[2])

  # Pre-compute some parts of the log density function.

  sigma.inv <- solve(sigma)
  log.det <- sum(log(eigen(sigma, only.values=TRUE)$values))
  leading.term <- -ndim/2 * log(2*pi) - log.det/2

  # Define log density and its gradient.

  log.density <- function(x) {
    stopifnot(length(x)==ndim)
    dist <- mahalanobis(as.vector(x), mean, sigma.inv, inverted=TRUE)
    return(leading.term-0.5*dist)
  }

  grad.log.density <- function(x) {
    stopifnot(length(x)==ndim)
    as.vector(-sigma.inv %*% array(x-mean, c(ndim,1)))
  }

  # overdispersed initial point

  initial <- function() {
    as.numeric(mvtnorm::rmvnorm(1, mean, 25/sqrt(ndim)*sigma))
  }

  mean.log.dens <- log.density(mean) - 0.5 * ndim

  # Define name and name.expression.

  if (is.null(rho)) {
    name <- sprintf('N%d', length(mean))
    name.expression <- sprintf('N[%d]', length(mean))
  } else {
    name <- sprintf('N%d,rho=%g', length(mean), rho)
    name.expression <- sprintf('N[%d](rho==%g)', length(mean), rho)
  }

  # Call make.dist with the functions we have defined.

  ds <- make.dist(ndim, name, name.expression, log.density,
                  grad.log.density, initial=initial, mean=mean, cov=sigma,
                  mean.log.dens=mean.log.dens)
}

# Define some Gaussians I use repeatedly.

N2weakcor.dist <- make.gaussian(c(0,0), rho=0.8)
N4poscor.dist <- make.gaussian(c(1,2,3,4), rho=0.999)
N4negcor.dist <- make.gaussian(c(1,2,3,4), rho=-0.3329)


#### RADFORD NEAL'S FUNNEL DISTRIBUTION ####

# See ?funnel.dist and R. M. Neal (2003), "Slice Sampling," Annals
# of Statistics 31(3):705-767, p. 732

funnel.dist <- make.dist(10, 'funnel', mean=rep(0,10),
  log.density = function(x) {
    dnorm(x[1],0,3,log=TRUE) + sum(dnorm(x[2:10],0,exp(0.5*x[1]), log=TRUE))
  },
  grad.log.density = function(x) {
    d <- -x/c(3^2,rep(exp(x[1]),9))
    d[1] <- d[1] + 0.5 * exp(-x[1]) * sum(x[2:10]^2) - 4.5
    return(d)
  }
)

#### GELMAN'S EIGHT SCHOOLS EXAMPLE ####

# See BDA, Gelman et al, sec 5.5, p. 138.

# This is a ten-dimensional multilevel model, with the first
# coordinate being a hyperparameter for the group means (mu), the
# second being a hyperparameter for the log of the variance of the
# group means (logtau), and the third through tenth (theta) being
# group means.  The schools.ybar and schools.var are means and
# variances of samples drawn from Gaussians with means of the
# corresponding thetas.

# This is expressed as a function that evaluates itself to keep the
# namespace clean.

schools.dist <- (function() {
  schools.ybar <- c(28, 8, -3, 7, -1, 1, 18, 12)
  schools.var <- c(15, 10, 16, 11, 9, 11, 10, 18)^2

  # Define log density and its gradient.

  log.density <- function(x) {
    stopifnot(all(is.finite(x)))
    mu <- x[1]
    logtau <- x[2]
    theta <- x[3:10]
    loglik <- sum( -0.5 * (schools.ybar-theta)^2/schools.var ) +
      sum( -0.5 * (theta-mu)^2/exp(2*logtau) ) - 8*logtau + logtau
    if (!is.finite(loglik))
      loglik <- -Inf
    return(loglik)
  }

  grad.log.density <- function(x) {
    mu <- x[1]
    logtau <- x[2]
    theta <- x[3:10]
    dtheta <- -(theta-schools.ybar)/schools.var -(theta-mu)*exp(-2*logtau)
    dmu <- -sum((mu-theta)/exp(2*logtau))
    dlogtau <- sum( (theta-mu)^2 ) * exp(-2*logtau) - 8 + 1
    return(c(dmu, dlogtau, dtheta))
  }

  # Define true means and covariances of parameters based on a long,
  # believed-to-be-reliable MCMC simulation.

  schools.mean <- c(7.831, 1.518, 11.406,  7.785, 5.708,
                   7.591, 5.071, 5.982, 10.683, 8.441)

  schools.cov <- rbind(
    c(27.84606,0.00152,20.12,15.9146,19.16,17.352,14.00,16.22,16.11,18.799),
    c(0.00152,1.10030,3.08,-0.0773,-2.06,-0.276,-2.14,-1.57,2.24,0.474),
    c(20.12353,3.08310,70.44,13.9797,7.83,11.875,4.87,7.45,22.21,14.313),
    c(15.91464,-0.07727,13.98,42.5085,15.00,11.642,11.39,11.67,10.75,12.722),
    c(19.16028,-2.05548,7.83,14.9955,69.78,14.533,16.74,17.99,9.09,14.027),
    c(17.35156,-0.27648,11.87,11.6419,14.53,47.133,12.05,13.18,11.14,16.373),
    c(14.00100,-2.13530,4.87,11.3909,16.74,12.045,39.25,14.64,5.04,8.935),
    c(16.22215,-1.57046,7.45,11.6657,17.99,13.176,14.64,48.55,8.34,12.196),
    c(16.11171,2.23727,22.21,10.7545,9.09,11.139,5.04,8.34,47.23,15.378),
    c(18.79947,0.47433,14.31,12.7219,14.03,16.373,8.93,12.20,15.38,62.346))

  # From a simulation of length 300k from double.sample.
  # 95% CI: (-16.8,-16.3)

  mean.log.dens <- -16.5

  # Call make.dist to glue the parts together.

  return(make.dist(10, '8-schools', 'plain("Eight Schools")',
                   log.density, grad.log.density,
                   mean=schools.mean, cov=schools.cov,
		   mean.log.dens=mean.log.dens))
})()


#### ROBERTS AND ROSENTHAL'S CONE DISTRIBUTION ####

# See ?make.cone.dist for more information.

make.cone.dist <- function(ndim) {
  make.c.dist(ndim, sprintf("cone%d", ndim), "cone_log_dens",
              name.expression=sprintf("cone[%d]", ndim),
              mean=rep(0,ndim))
}


#### RANDOM MULTIMODAL DISTRIBUTIONS ####

# See ?make.multimodal.dist for more information.

make.multimodal.dist <- function(nmodes, ndim, cube.size) {

  # Set random seed to fixed value.

  runif(1)  # force .Random.seed to exist
  seed <- .Random.seed
  set.seed(17)

  # Draw random modes and compute the overall mean.

  modes <- array(runif(ndim*nmodes,0,cube.size),c(nmodes,ndim))
  mu <- as.numeric(colMeans(modes))

  # Define log density and its gradient.

  L <- function(x) {
    if (is.matrix(x))
      x <- as.vector(x)
    stopifnot(is.vector(x) && length(x)==ndim)
    y <- logsumexp(-0.5 * rowSums(sweep(modes, 2, x)^2))
    return(y)
  }

  dL <- function(x) {
    if (is.matrix(x))
      x <- as.vector(x)
    stopifnot(is.vector(x) && length(x)==ndim)
    colSums(exp(-0.5*rowSums(sweep(modes, 2, x)^2)) * sweep(modes, 2, x)) /
      exp(L(x))
  }

  # Restore random seed.

  .Random.seed <<- seed

  # Call make.dist to define the dist object, stash the modes for
  # convenience, and return the object.

  name <- sprintf("mixture%d*%d", ndim, nmodes)
  ds <- make.dist(ndim, name, log.density=L, grad.log.density=dL, mean=mu)
  ds$modes <- modes
  return(ds)
}


#### UNCORRELATED GAMMAS ####

# shape is a vector of shapes, scale is a vector of scales.  They must
# be the same length.  See ?make.mv.gamma.dist for more information.
# ?make.dist has a simpler example of Gamma distribution objects.

make.mv.gamma.dist <- function(shape, scale=rep(1, length(shape))) {
  stopifnot(length(shape)==length(scale))
  name <- sprintf("Gamma%d", length(shape))
  name.expression <- sprintf("Gamma[%d]", length(shape))

  # Precompute the normalizing constant, then define the distribution
  # object.

  Z <- lgamma(shape) + shape*log(scale)
  ds <- make.dist(length(shape), name, name.expression,
    log.density=function(x)
      ifelse(any(x<0), -Inf, sum((shape-1)*log(x) - x/scale - Z)),
    grad.log.density=function(x) ifelse(x<0, Inf, (shape-1)/x - 1/scale),
    mean=shape*scale,
    cov=diag(shape*scale^2, nrow=length(shape)))

  # Save the parameters for convenience.

  ds$shape <- shape
  ds$scale <- scale

  # Undocumented API for specifying lower bound on support to
  # methods that know about it.

  ds$lower.bound <- rep(0, length(shape))

  return(ds)
}
