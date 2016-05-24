summary.Bayesianbetareg <-
function(object, ...){
  
  se <- object$desv
  
  TAB <- cbind( Coefficient = coef(object),
                Desv. = se,
                L.CredIntv = object$interv[,1],
                U.CredIntv = object$interv[,2]
  )
  
  colnames(TAB) <- c("Estimate", "Est.Error", "L.CredIntv",  "U.CredIntv")
  
  beta.resid <- betaresiduals(object$Y,object$X,object)
  criteria <- criteria(object$X,beta.resid)
  
  res <- list(call=object$call, coefficients=TAB, Deviance=criteria$Deviance, AIC=criteria$AIC, BIC=criteria$BIC)  
  
  class(res) <- "summary.Bayesianbetareg"
  res  
}
