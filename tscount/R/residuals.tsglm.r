residuals.tsglm <- function(object, type=c("response", "pearson", "anscombe"), ...){
  type <- match.arg(type)
  if(type=="pearson"){
    #Standard deviation of the conditional distribution (cf. marcal.tsglm):
    if(object$distr=="poisson") sddistr <- function(meanvalue, distrcoefs) sqrt(meanvalue)
    if(object$distr=="nbinom") sddistr <- function(meanvalue, distrcoefs) sqrt(meanvalue + meanvalue^2/distrcoefs)
    result <- object$residuals/sddistr(meanvalue=fitted(object), distrcoefs=object$distrcoefs)
  }
  if(type=="anscombe"){
    y <- object$response
    mu <- fitted(object)
    if(object$distr=="poisson") result <- 3/2*(y^(2/3)-mu^(2/3))/mu^(1/6)
    if(object$distr=="nbinom") result <- (3*object$distrcoefs*((1+y/object$distrcoefs)^(2/3) - (1+mu/object$distrcoefs)^(2/3)) + 3*(y^(2/3)-mu^(2/3))) / (2*(mu^2/object$distrcoefs + mu)^(1/6))
  }else{
    result <- object$residuals
  }
  return(result)
}
