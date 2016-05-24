predict.prm <-
function(object, newdata, ...){
  
# PREDICT.PRM predicts response based on a Partial Robust M regression model 
#     Inputs: object, a "prm" class Partial Robust M regression object 
#             newdata, an optional data matrix or frame containing a new set of cases 
#
# Written by Sven Serneels, BASF Corp, December 2013. 
 
  if(!(class(object)=="prm")){stop("The PRM predict function only applies to prm class objects")}
  b <- coef(object)
  if(missing(newdata)){fitted.values <- object$fitted.values}
  else{
    if(ncol(newdata)==length(b)){Xn <- newdata} 
    else{
      formula <- object$input$formula
      mt <- terms(formula, data=object$input$X0)
      ic <- attr(mt, "intercept")
      if (ic==0){
        Xn <- model.matrix(mt, newdata)
      } else{
        Xn <- model.matrix(mt, newdata)[,-1]
      } 
    }  
	if (is.vector(Xn)==TRUE){
		Xn <- t(as.matrix(Xn))
	}
    fitted.values <- (as.matrix(Xn)%*%b + object$intercept)[,1]
  }
  return(fitted.values)
}
