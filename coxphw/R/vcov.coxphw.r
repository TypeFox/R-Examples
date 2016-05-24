vcov.coxphw <- function(object, ...)
{
  if (!is.null(object$betafix)) { 
    cat(paste("The following variables were not estimated in this model",               
              paste(names(object$coefficients)[!is.na(object$betafix)], collapse=", "), 
              "\n\n"), sep="") }    
  return(object$var)
}