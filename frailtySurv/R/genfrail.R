genfrail <- function(# Number of clusters and cluster sizes
                     N = 300,
                     # If K is int, then fixed cluster sizes
                     # If K is numeric vector, then sizes provided
                     # Otherwise draw random
                     # can be one of: "poisson", "pareto", "uniform"
                     K = 2, 
                     # cluster size distributian params:
                     # K-truncated Poisson: c(lambda, truncated value)
                     # Discrete Pareto: c(alpha, max value inclusive)
                     # Uniform: c(lower, upper) inclusive
                     K.param = c(2, 0), 
                     
                     # Covariate coefficients
                     beta = c(log(2)),
                     
                     # Frailty distribution and parameter vector
                     # Can be one of: "gamma", "pvf", "lognormal", "invgauss", "posstab", "none"
                     frailty = "gamma",
                     # Frailty distribution parameter vector
                     theta = c(2), 
                     
                     # Covariate distribution and params
                     # covar.distr can be one of: "normal", "uniform", "zero"
                     covar.distr = "normal", 
                     covar.param = c(0,1),
                     covar.matrix = NULL,
                     
                     # Censoring distribution and parameters vector
                     # censor.distr can be one of: "normal", "lognormal", "uniform", "none"
                     censor.distr = "normal",
                     censor.param = c(130,15),
                     censor.rate = NULL, # If specified, overrides censor.mu
                     
                     # Only one of these needs to be specified
                     # Order of preference is: Lambda_0_inv, Lambda_0, lambda_0
                     lambda_0 = NULL, #i.e. function(t, tau=4.6, C=0.01) (tau*(C*t)^tau)/t,
                     Lambda_0 = NULL, #i.e. function(t, tau=4.6, C=0.01) (C*t)^tau,
                     Lambda_0_inv = NULL, #i.e. function(t, tau=4.6, C=0.01) (t^(1/tau))/C,
                     
                     # Round time to nearest round.base
                     round.base = NULL
) {
  Call <- match.call()
  
  # Determine cluster sizes
  if (is.numeric(K) && length(K) == 1) {
    cluster.sizes <- c(rep(K, N));
  } else if (is.numeric(K) && length(K) > 1) {
    if (length(K) != N)
      stop("length(K) must equal N when providing cluster sizes")
    cluster.sizes <- K
  } else if (K == "poisson") {
    cluster.sizes <- rtpois(N, K.param[1], K.param[2])
  } else if (K == "pareto") {
    cluster.sizes <- rtzeta(N, K.param[1], K.param[2], K.param[3])
  } else if (K == "uniform") {
    cluster.sizes <- round(runif(N, K.param[1], K.param[2]))
  } else {
    stop("Wrong value for K. Must be int, string, or numeric vector")
  }
  
  NK <- sum(cluster.sizes)
  
  # Generate the covariates
  p <- length(beta)
  Z <- matrix(0, nrow=NK, ncol=p)
  for (j in 1:p) {
    if (!is.null(covar.matrix)) {
      Z <- covar.matrix
    } else if (covar.distr == "normal") {
      Z[, j] <- rnorm(NK, covar.param[1], covar.param[2])
    } else if (covar.distr == "uniform") {
      Z[, j] <- runif(NK, covar.param[1], covar.param[2])
    } else if (covar.distr == "zero") {
      Z[, j] <- rep(0, NK)
      covar.param <- NULL
    }
  }
  
  # Generate the frailty for each cluster
  if (frailty %in% names(rfrailty)) {
    if (!all((theta > lb.hard.frailty[[frailty]])&(theta < ub.hard.frailty[[frailty]])))
        stop(frailty, " frailty distribution parameters must be in the range (",
             lb.hard.frailty[[frailty]], ",", ub.hard.frailty[[frailty]], ")")
    cluster.frailty <- rfrailty[[frailty]](N, theta)
  } else if (frailty == "none") {
    cluster.frailty <- rep(1, N)
  } else {
    stop("Wrong value for frailty: ", frailty)
  }
  
  w = rep(cluster.frailty, cluster.sizes)
  rate <- w*exp(Z%*%beta)
  
  if (sum(c(!is.null(lambda_0), !is.null(Lambda_0), !is.null(Lambda_0_inv)) > 1))
    stop("Must specify only one of: lambda_0, Lambda_0, or Lambda_0_inv")
  
  # Generate uniform random samples for the survival times
  u <- runif(NK)
  
  # With the inverse cumulative hazard, apply Bender's method
  if (!is.null(Lambda_0_inv)) {
    fail.time <- Lambda_0_inv(-log(u)/rate)
    
    # hazard attribute, for summary
    hazard <- list(Lambda_0_inv=Lambda_0_inv)
    
    # Inverse to get the original CBH
    Lambda_0 <- Vectorize(function(y) uniroot(function(x) 
      Lambda_0_inv(x) - y, lower=0, upper=100, extendInt="yes")$root)
  # Otherwise, use the method described by Crowther
  } else {
    
    if (is.null(Lambda_0) && !is.null(lambda_0)) {
      # numeric integral nested in the root finding
      Lambda_0 <- Vectorize(function(t) {
        if (t <= 0) {
          return(0)
        }
        integrate(lambda_0, 0, t, subdivisions = 1000L)$value
      })
      hazard <- list(lambda_0=lambda_0)
    } else if (!is.null(Lambda_0)) {
      hazard <- list(Lambda_0=Lambda_0)
    } else {
      warning("Using a default baseline hazard. Did you forget to pass Lambda_0?")
      Lambda_0 <- function(t, tau=4.6, C=0.01) (C*t)^tau
      hazard <- list(Lambda_0=Lambda_0)
    }
    
    # TODO: first try to invert the function analytically, maybe using Ryacas?
    
    # if inversion fails, then need to find root for each t
    # Similar to the method described by Crowther, 
    # except solve: log(S(t)) - log(u) = 0
    fail.time <- sapply(1:(NK), function(ij) {
      fn <- function(t) {
        (-Lambda_0(t)*rate[ij]) - log(u[ij])
      }
      # nleqslv(100, fn, global="dbldog", control=list(allowSingular=TRUE))$x
      uniroot(fn, lower=0, upper=100, extendInt="yes")$root
    })
  }
  
  if (censor.distr == "normal") {
    censor.density <- dnorm
    censor.random <- rnorm
    censor.quantile <- qnorm
    censor.mu <- censor.param[1]
    censor.sigma <- censor.param[2]
  } else if (censor.distr == "lognormal") {
    censor.density <- dlnorm
    censor.random <- rlnorm
    censor.quantile <- qlnorm
    censor.sigma <- sqrt(log(1 + censor.param[2]^2/censor.param[1]^2))
    censor.mu <- log(censor.param[1]) - censor.sigma^2/2
  } else if (censor.distr == "uniform") {
    censor.density <- dunif
    censor.random <- runif
    censor.quantile <- qunif
    censor.lower <- censor.param[1]
    censor.upper <- censor.param[2]
    stopifnot(censor.lower < censor.upper)
  } else if (censor.distr == "none") {
    warning("Censoring distribution is set to none")
  } else {
    stop("Censoring distribution must be normal or lognormal")
  }
  
  if (censor.distr == "none") {
    obs.status <- rep(1, NK)
    obs.time <- fail.time
  } else {
    if (!is.null(censor.rate)) {
      stopifnot(0 <= censor.rate && censor.rate <= 1)
      
      # Empirical survivor and hazard functions
      esurv <- ecdf(fail.time)
      efail <- function(t) 1 - esurv(t)
      
      if (censor.distr == "uniform") {
        interval <- censor.upper - censor.lower
        # Split the integral up into 2 parts since it may actually be ~0 for 
        # most of the interval. Make this split at the 50th percentile
        mu <- uniroot(function(mu) {
            lower <- mu - interval/2
            upper <- mu + interval/2
            censor.rate - 
              (integrate(function(t) efail(t)*censor.density(t, lower, upper), 
                         -Inf, censor.quantile(0.5, lower, upper), subdivisions=1e3)$value +
                 integrate(function(t) efail(t)*censor.density(t, lower, upper), 
                           censor.quantile(0.5, lower, upper), Inf, subdivisions=1e3)$value)
          }, lower=0, upper=100, extendInt="up")$root
        censor.lower <- mu - interval/2
        censor.upper <- mu + interval/2
      } else {
        # Split the integral up into 2 parts since it may actually be ~0 for 
        # most of the interval. Make this split at the 50th percentile
        censor.mu <- uniroot(function(mu) 
          censor.rate - 
            (integrate(function(t) efail(t)*censor.density(t, mu, censor.sigma), 
                       -Inf, censor.quantile(0.5, mu, censor.sigma), subdivisions=1e3)$value +
               integrate(function(t) efail(t)*censor.density(t, mu, censor.sigma), 
                         censor.quantile(0.5, mu, censor.sigma), Inf, subdivisions=1e3)$value),
          lower=0, upper=100, extendInt="up")$root
      }
    }
    
    # Non-informative right censoring
    if (censor.distr == "uniform") {
      censor.time <- censor.random(NK, censor.lower, censor.upper)
    } else {
      censor.time <- censor.random(NK, censor.mu, censor.sigma)
    }
    obs.status <- sign(fail.time <= censor.time)
    obs.time <- pmin(fail.time, censor.time)
  }
  
  if (is.numeric(round.base)) {
    obs.time <- round.base*floor(obs.time/round.base + 0.5)
  }
  
  # Covariates are named Z1, Z2, ...
  colnames(Z) <- paste("Z", 1:p, sep="")
  
  dat <- data.frame(family=rep(1:N, cluster.sizes),
                    rep=unlist(lapply(cluster.sizes, function(m) rep(1:m))), 
                    time=obs.time, 
                    status=obs.status,
                    Z)
  class(dat) <- append("genfrail", class(dat))
  
  # These are mostly needed for printing the summary and simulations
  attributes(dat) <- append(attributes(dat), list(
    created=Sys.time(),
    beta=setNames(beta, paste("beta.", 1:length(beta), sep="")),
    theta=setNames(theta, paste("theta.", 1:length(theta), sep="")),
    frailty=frailty,
    covar.distr=covar.distr,
    covar.param=covar.param,
    Lambda_0=Lambda_0,
    hazard=hazard,
    call=Call
  ))
  
  dat
}
