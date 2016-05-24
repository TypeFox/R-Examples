#'@export
predict.DStree <- function(object,data,...){
  
  fitr<-object
  class(fitr)<-"rpart"
  pred <- predict(fitr,data,type="matrix")
  col<-ncol(pred)/2
  return(list(MedSurv=predict(fitr,data,type="vector"),
          Haz=pred[,1:col],Surv=pred[,(col+1):(2*col)]))
  
}