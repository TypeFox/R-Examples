predict.sprm <-
function(object,newdata, ...){
  
# PREDICT.SPRM predicts response based on a Sparse Partial Robust M regression model 
#     Inputs: object, a "sprm" class Partial Robust M regression object 
#             newdata, an optional data matrix or frame containing a new set of cases 
#
# Written by Sven Serneels, BASF Corp, January 2014. 

  if(!(class(object)=="sprm")){stop("The SPRM predict function only applies to sprm class objects")}
  b <- coef(object)
  Xn <- NULL
  if(missing(newdata)){fitted.values <- object$fitted.values}
  else{
    if(ncol(newdata)==length(b)){Xn <- newdata} 
    else{
      formula <- object$input$formula
      mf <- formula(paste(as.character(formula)[c(1,3)], sep=""))
      mt <- terms(mf, data=object$input$X0)      
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
    fitted.values <- as.vector(as.matrix(Xn)%*%b + object$intercept)
  }
  return(fitted.values)
}
