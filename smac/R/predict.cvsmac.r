predict.cvsmac <- function(object,new.x = NULL,...){
  if (is.null(new.x)) {new.x=object$model$x}
  
  # calling predict.smac
    
  best.lambda <- object$best.lambda
  b <- predict.smac(object$model, new.x = new.x, lambda = best.lambda)
  best <- list(new.x=new.x, best.lambda = best.lambda, best.beta0 = b$fitted.beta0, best.beta = b$fitted.beta, best.pred.y = b$pred.y, best.pred.prob = b$pred.prob)
  
  # Returning values  
  this.call = match.call()
  best$call <- this.call
  
  return(best)
}
