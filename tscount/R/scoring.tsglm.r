scoring <- function(object, ...) UseMethod("scoring")

scoring.tsglm <- function(object, cutoff=1000, ...){
  #Density functton of the conditional distribution:
  if(object$distr=="poisson") ddistr <- function(x, meanvalue, distrcoefs) dpois(x, lambda=meanvalue)
  if(object$distr=="nbinom") ddistr <- function(x, meanvalue, distrcoefs) dnbinom(x, mu=meanvalue, size=distrcoefs[["size"]])  
  #Cumulative distribution function of the conditional distribution (cf. marcal.tsglm):
  if(object$distr=="poisson") pdistr <- function(q, meanvalue, distrcoefs) ppois(q, lambda=meanvalue)
  if(object$distr=="nbinom") pdistr <- function(q, meanvalue, distrcoefs) pnbinom(q, mu=meanvalue, size=distrcoefs[["size"]])
  #Standard deviation of the conditional distribution (cf. marcal.tsglm):
  if(object$distr=="poisson") sddistr <- function(meanvalue, distrcoefs) sqrt(meanvalue)
  if(object$distr=="nbinom") sddistr <- function(meanvalue, distrcoefs) sqrt(meanvalue + meanvalue^2/distrcoefs[["size"]])
  n <- object$n_eff 
  logarithmic <- quadratic <- spherical <- rankprob <- dawseb <- normsq <- sqerror <- numeric(n) #scores
  for(t in 1:n){
    #auxiliary objects, which are overwritten in each step:
      y <- object$response[t]
      mu <- fitted(object)[t]
      sigma <- sddistr(meanvalue=mu, distrcoefs=object$distrcoefs)
      p_y <- ddistr(y, meanvalue=mu, distrcoefs=object$distrcoefs)
      quadrat_p <- sum(ddistr(0:cutoff, meanvalue=mu, distrcoefs=object$distrcoefs)^2)
    #computation of the scores:
      logarithmic[t] <- - log(p_y)
      quadratic[t] <- - 2*p_y + quadrat_p
      spherical[t] <- - p_y/sqrt(quadrat_p)
      rankprob[t] <- sum((pdistr(0:cutoff, meanvalue=mu, distrcoefs=object$distrcoefs) - as.integer(y <= 0:cutoff))^2)
      sqerror[t] <- (y-mu)^2
      normsq[t] <- sqerror[t]/sigma^2 
      dawseb[t] <- normsq[t] + 2*log(sigma)
  }
  result <- c(
    logarithmic=mean(logarithmic),
    quadratic=mean(quadratic),
    spherical=mean(spherical),
    rankprob=mean(rankprob),
    dawseb=mean(dawseb),
    normsq=mean(normsq),
    sqerror=mean(sqerror)
  )
  return(result)
}
