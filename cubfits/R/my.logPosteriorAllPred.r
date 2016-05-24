### This either returns a log posterior vector or sum of the vector.
###
### These functions are for all predicting genes.

### Get the specific function according to the options.
get.my.logPosteriorAllPred <- function(model.Phi){
  if(!any(model.Phi[1] %in% .CF.CT$model.Phi)){
    stop("model.Phi is not found.")
  }
  ret <- eval(parse(text = paste("my.logPosteriorAllPred.",
                                 model.Phi[1], sep = "")))
  assign("my.logPosteriorAllPred", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.logPosteriorAllPred().


### Function to calculate complete log-posterior for
### (phi, b, mu.Phi, sigma.Phi.sq) given y and n
my.logPosteriorAllPred.lognormal <- function(phi, y, n, b, p.Curr,
    reu13.df = NULL){
  ### Dispatch.
  nu.Phi <- p.Curr[1]
  sigma.Phi <- p.Curr[2]

  ### Return.
  ret <- .cubfitsEnv$my.logdmultinomCodAllR(b, phi, y, n, reu13.df = reu13.df) +
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
} # End of my.logPosteriorAllPred.lognormal().

### Function for log mixture normal.
my.logPosteriorAllPred.logmixture <- function(phi, y, n, b, p.Curr,
    reu13.df = NULL){
  ### Dispatch.
  paramlog <- p.Curr

  ### Return
  ret <- .cubfitsEnv$my.logdmultinomCodAllR(b, phi, y, n, reu13.df = reu13.df) +
         dlmixnorm(log(phi), paramlog, log = TRUE)
  ret
} # End of my.logPosteriorAllPred.logmixture().
