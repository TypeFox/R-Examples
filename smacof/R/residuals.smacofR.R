#residuals of distances for smacof objects

residuals.smacofR <- function(object, ...)
{
  # object of class smacof
  return(as.matrix(object$obsdiss - object$confdiss))
}
