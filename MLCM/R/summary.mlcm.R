`summary.mlcm` <-
function(object, 
	digits = max(3, getOption("digits") - 4), ...){
  z <- object
  ans <- list()
  ans$pscale <- z$pscale	
  ans$se <- if(z$method == "glm") 
  		coef(summary(z$obj))[, 2] else
  		object$se  
  ans$sigma <- z$sigma
  ans$logLik <- logLik(z)[1]
  ans$aic <- AIC(z)	
  ans$link <- z$link
  ans$method <- z$method
  ans$par <- z$par
  ans$model <- z$model
  ans$formula <- z$formula
  class(ans) <- "summary.mlcm"
  ans
}

