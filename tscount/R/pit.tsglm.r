pit <- function(object, ...) UseMethod("pit")

pit.tsglm <- function(object, bins=10, ...){
  #Cumulative distribution function of the conditional distribution (cf. marcal.tsglm):
  if(object$distr=="poisson") pdistr <- function(q, meanvalue, distrcoefs) ppois(q, lambda=meanvalue)
  if(object$distr=="nbinom") pdistr <- function(q, meanvalue, distrcoefs) pnbinom(q, mu=meanvalue, size=distrcoefs[["size"]])
  n <- object$n_eff
  u <- seq(0, 1, length=bins+1)
  pit <- numeric(length(u))
  for(i in 1:n){
    P_x <- pdistr(object$response[i], meanvalue=fitted(object)[i], distrcoefs=object$distrcoefs)
    if(object$response[i]!=0){
      P_x_1 <- pdistr(object$response[i]-1, meanvalue=fitted(object)[i], distrcoefs=object$distrcoefs)
    }else{
      P_x_1 <- 0
    }
    pit <- pit + punif(u, P_x_1, P_x)/n
  }
  histo <- list(breaks=u, counts=diff(pit)*n, density=diff(pit)*bins, mids=(u[-(bins+1)]+u[-1])/2, xname="PIT", equidits=TRUE)
  class(histo) <- "histogram"
  #simconfint <- if(ci>0 && ci<1) (n/bins+c(-1,+1)*qnorm(1-(1-ci)/bins/2)*sqrt(n*(1/bins)*(1-1/bins)))/(n/bins) else NULL #simultaneous confidence band of level ci (normal approximation) for the histogram bars under the assumption of iid U(0,1) PIT values 
  plot_args <- modifyList(list(main="Non-randomized PIT histogram", xlab="Probability integral transform", ylab="Density", freq=FALSE, ylim=range(0, histo$density)), list(...)) #the default arguments can be overriden by those provided in the ... argument
  do.call("plot", args=c(list(x=histo), plot_args))
  #if(ci>0 && ci<1) abline(h=simconfint, lty="dashed", col=ci.col)
}
