summary.Gammareg <-
function(object, ...){
  

  
  TAB <- cbind( Coefficient = coef(object),
                L.CredIntv = object$interv[,1],
                U.CredIntv = object$interv[,2]
  )
  
  colnames(TAB) <- c("Estimate", "L.Intv",  "U.Intv")
  
  
  res <- list(call=object$call, coefficients=TAB, covB=object$desvB, covG=object$desvG, 
  AIC=object$AIC, iteration=object$iteration, convergence=object$convergence)  
  
  class(res) <- "summary.Gammareg"
  res  
}
