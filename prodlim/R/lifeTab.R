#  These functions extract the number of subjects atRisk and the number of
#  events at given times from the object and binds it together with 
#  quantities like survival prob, cuminc, standard errors, etc. which can 
#  simply be evaluated at the requested times.

lifeTab <- function(object,...){
  if(NROW(object$model.response)<=0) stop("No response found") #  to avoid seg faults
  dummy <- 1
  class(dummy) <- object$model
  UseMethod("lifeTab",object=dummy)
}






  
