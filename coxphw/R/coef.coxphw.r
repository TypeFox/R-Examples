coef.coxphw <- function(object, ...)
{
  #if (object$template=="PH")  { coef <- coef(object) } else
  if (object$template=="PH")  { coef <- object$coefficients } else
  if (object$template!="PH") { 
    coef <-as.vector(t(object$coefficients))
    names(coef) <- names(object$coefficients)
  }   
  
  return(coef)
}