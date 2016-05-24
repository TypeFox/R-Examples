autojags <- function(object, n.iter=1000, n.thin=1, Rhat=1.1, n.update=2, refresh=n.iter/50, 
    progress.bar = "text",...)
{
  n.burnin = object$n.iter
  n.thin.auto <- max( 1, floor( ( n.iter - n.burnin )/1000 ) )
  n.thin <- ifelse(n.thin > n.thin.auto, n.thin, n.thin.auto)

  if(any(!class(object) %in% c("rjags","rjags.parallel"))) stop("model must be a rjags object")
    object <- update(object, n.iter=n.iter, n.thin=n.thin, 
                      refresh=refresh, progress.bar = progress.bar,...)
    check <- any(object$BUGSoutput$summary[,"Rhat"] > Rhat)
    if (check){
      count <- 1
      while (check & (count < n.update)) {
          object <- update(object, n.iter=n.iter, n.thin=n.thin, 
                      refresh=refresh, progress.bar = progress.bar, ...)
          count <- count + 1
          check <- any(object$BUGSoutput$summary[,"Rhat"] > Rhat)
      }
    }
    return(object)
}
