`logLik.mws` <-
function(object,...)
#returns the log-likelihood value of an object of class mws
{
  llik <- object$likelihood[length(object$likelihood)]
  cat("Log-Likelihood:", llik,"\n")
}

