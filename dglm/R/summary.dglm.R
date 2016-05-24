summary.dglm <- function(object, dispersion=NULL, correlation = FALSE, ...)
{
  #  Summarize a double glm
  #  GKS  7 Jan 98
  sm <- stats::summary.glm(object,dispersion = dispersion, correlation = correlation, ...)
  sd <- stats::summary.glm(object$dispersion.fit,dispersion = 2, correlation = correlation, ...)
  ans <- structure(c(sm, list(dispersion.summary = sd, outer.iter = object$iter,
                              m2loglik = object$m2loglik)), 
                   class = c("summary.dglm","summary.glm"))
  
  class(ans) <- "summary.dglm"
  return(ans)
}
