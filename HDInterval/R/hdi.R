# This file has the S3 generic 'hpdi' function and a series of methods.

hdi <- function(object, credMass=0.95, ...) UseMethod("hdi")

hdi.default <- function(object, credMass=0.95, ...) {
  if(!is.numeric(object))
    stop(paste("No applicable method for class", class(object)))
  checkCredMass(credMass)
  result <- hdiVector(object, credMass=credMass)
  # names(result) <- c("lower", "upper")
  attr(result, "credMass") <- credMass
  return(result)
}

hdi.matrix <- function(object, credMass=0.95, ...) {
  checkCredMass(credMass)
  result <- apply(object, 2, hdiVector, credMass=credMass)
  attr(result, "credMass") <- credMass
  return(result)
}

hdi.data.frame <- function(object, credMass=0.95, ...) {
  checkCredMass(credMass)
  result <- sapply(object, hdiVector, credMass = credMass)
  attr(result, "credMass") <- credMass
  return(result)
}

hdi.mcmc.list <- function(object, credMass=0.95, ...)
  hdi.matrix(as.matrix(object), credMass=credMass, ...)

hdi.mcmc <- function(object, credMass=0.95, ...)
  hdi.matrix(as.matrix(object), credMass=credMass, ...)

hdi.bugs <- function(object, credMass=0.95, ...)
  hdi.matrix(object$sims.matrix, credMass=credMass, ...)

hdi.rjags <- function(object, credMass=0.95, ...)
  hdi.matrix(object$BUGSoutput$sims.matrix, credMass=credMass, ...)

hdi.runjags <- function(object, credMass=0.95, ...)
  hdi.mcmc.list(object$mcmc, credMass=credMass, ...)

hdi.jagsUI <- function(object, credMass=0.95, ...)
  hdi.mcmc.list(object$samples, credMass=credMass, ...)




