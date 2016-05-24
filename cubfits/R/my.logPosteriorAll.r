### This either returns a log posterior vector or sum of the vector.
###
### These functions are for all genes.

### Get the specific function according to the options.
get.my.logPosteriorAll <- function(model.Phi){
  if(!any(model.Phi[1] %in% .CF.CT$model.Phi)){
    stop("model.Phi is not found.")
  }
  if(!.CF.CONF$estimate.bias.Phi){
    ret <- eval(parse(text = paste("my.logPosteriorAll.",
                                   model.Phi[1], sep = "")))
  } else{
    ret <- eval(parse(text = paste("my.logPosteriorAll.",
                                   model.Phi[1], "_bias", sep = "")))
  }
  assign("my.logPosteriorAll", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.logPosteriorAll().


### Function to calculate complete log-posterior for
### (phi, b, sigmaWsq, mu.Phi, sigma.Phi.sq) given y, n, and phi.Obs
my.logPosteriorAll.lognormal <- function(phi, phi.Obs, y, n, b, p.Curr,
    reu13.df = NULL){
  ### Dispatch.
  sigmaW <- p.Curr[1]
  nu.Phi <- p.Curr[2]
  sigma.Phi <- p.Curr[3]

  ### Return.
  ret <- dlnorm(phi.Obs, log(phi), sigmaW, log = TRUE) +
         .cubfitsEnv$my.logdmultinomCodAllR(b, phi, y, n, reu13.df = reu13.df) +
         dlnorm(phi, nu.Phi, sigma.Phi, log = TRUE)

  ### Note that prior for all Phi is lognormal(mu.Phi, sigma.Phi), and assume
  ### proper priors for mu.Phi (normal) and sigma.Phi (inverse gamma.)
  ### It is necessary to add this for proportional posterior probability.
  ### No need to add this for acceptance ratio in our random walk, since
  ### the acceptance ratio cancels out this term, and we assume non-informative
  ### prior for mu.Phi and sigma.Phi.sq.
  ### i.e. p(mu.Phi, sigma.Phi.sq) \propto 1/sigma.Phi.sq.
  ### See Gelman et al. (2003), p.75 for details.

  # ret <- ret +
  #        dnorm(mu.Phi, some.mean, some.sd, log = TRUE) +
  #        dgamma(1/sigma.Phi.sq, some.alpha, some.beta, log = TRUE) + sigma.Phi.sq^2

  ret
} # End of my.logPosteriorAll.lognormal().

### Function for log mixture normal.
my.logPosteriorAll.logmixture <- function(phi, phi.Obs, y, n, b, p.Curr,
    reu13.df = NULL){
  ### Dispatch.
  sigmaW <- p.Curr[1]
  paramlog <- p.Curr[-1]

  ### Return
  log.phi <- log(phi)
  ret <- dlnorm(phi.Obs, log.phi, sigmaW, log = TRUE) +
         .cubfitsEnv$my.logdmultinomCodAllR(b, phi, y, n, reu13.df = reu13.df) +
         dlmixnorm(log.phi, paramlog, log = TRUE)
  ret
} # End of my.logPosteriorAll.logmixture().


### Bias of Phi version.

### Function to calculate complete log-posterior for
### (phi, b, sigmaWsq, mu.Phi, sigma.Phi.sq) given y, n, and phi.Obs
my.logPosteriorAll.lognormal_bias <- function(phi, phi.Obs, y, n, b, p.Curr,
    reu13.df = NULL){
  ### Dispatch.
  sigmaW <- p.Curr[1]
  nu.Phi <- p.Curr[2]
  sigma.Phi <- p.Curr[3]
  bias.Phi <- p.Curr[4]

  ### Return.
  log.phi <- log(phi) + bias.Phi
  ret <- dlnorm(phi.Obs, log.phi, sigmaW, log = TRUE) +
         .cubfitsEnv$my.logdmultinomCodAllR(b, phi, y, n, reu13.df = reu13.df) +
         dlnorm(phi, nu.Phi, sigma.Phi, log = TRUE)
  ret
} # End of my.logPosteriorAll.lognormal_bias().

### Function for log mixture normal.
my.logPosteriorAll.logmixture_bias <- function(phi, phi.Obs, y, n, b, p.Curr,
    reu13.df = NULL){
  ### Dispatch.
  sigmaW <- p.Curr[1]
  paramlog <- p.Curr[-c(1, length(p.Curr))]
  bias.Phi <- p.Curr[length(p.Curr)]

  ### Return
  log.phi <- log(phi) + bias.Phi
  ret <- dlnorm(phi.Obs, log.phi, sigmaW, log = TRUE) +
         .cubfitsEnv$my.logdmultinomCodAllR(b, phi, y, n, reu13.df = reu13.df) +
         dlmixnorm(log.phi, paramlog, log = TRUE)
  ret
} # End of my.logPosteriorAll.logmixture_bias().
