marcal <- function(object, ...) UseMethod("marcal")

marcal.tsglm <- function(object, plot=TRUE, ...){
  #Cumulative distribution function of the conditional distribution:
  if(object$distr=="poisson") pdistr <- function(q, meanvalue, distrcoefs) ppois(q, lambda=meanvalue)
  if(object$distr=="nbinom") pdistr <- function(q, meanvalue, distrcoefs) pnbinom(q, mu=meanvalue, size=distrcoefs)
  xvalues <- min(object$ts):max(object$ts) #range of values could be extended in the future
  p_bar <- g_hat <- numeric(length(xvalues))
  n <- object$n_eff 
  for(i in seq(along=xvalues)){
    for(j in 1:n){
      p_bar[i] <- p_bar[i] + pdistr(xvalues[i], meanvalue=fitted(object)[j], distrcoefs=object$distrcoefs)/n
    }
    g_hat[i] <- sum(object$response<=xvalues[i])/n
  }
  result <- list(x=xvalues, y=p_bar-g_hat)
  if(plot){
    plot_args <- modifyList(list(main="Marginal calibration plot", xlab="Threshold value", ylab="Diff. of pred. and emp. c.d.f"), list(...)) #the default arguments can be overriden by those provided in the ... argument
    do.call("plot", args=c(list(result$x, result$y, type="l"), plot_args))
    abline(h=0, lty=3)
    invisible(result)
  }else{
    return(result)
  }
}
