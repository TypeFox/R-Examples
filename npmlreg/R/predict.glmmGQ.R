"predict.glmmGQ" <-
function(object, newdata,  type="link", ...){
  if (missing(newdata)){#Emp. Bayes Pred. (Aitkin, 96)
      if (type=="link") {
          return(round(object$ebp,digits=4))
      } else {
          return(round(object$fitted.values,digits=4))
      }
  } else {      # Analytical mean of compounded model, see Aitkin, Francis, Hinde 2005, p. 459.
      if (object$family$link!="log"){
          stop("Prediction for objetcs of class glmmGQ only supported for log link")
      }
      Terms<-  delete.response(terms(object$formula))
      X<-model.matrix(Terms, model.frame(Terms,newdata))
      pred<-    as.vector(X%*%matrix(object$coef[1:dim(X)[2]])+1/2* (object$coef[length(object$coef)])^2)
      if (type=="link"){
          rpred<-as.vector(pred)
      } else {
          rpred<- exp(as.vector(pred))
      }
      names(rpred)<-dimnames(X)[[1]]
      return(round(rpred,digits=4))
  }
}
