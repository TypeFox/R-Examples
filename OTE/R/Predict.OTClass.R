Predict.OTClass <-
function(Opt.Trees,XTesting, YTesting=NULL){
  
  if(is.null(YTesting)){
    Predictions <- predict(Opt.Trees$t.object, XTesting)
    mylist <- list("Predicted.Class.Labels"=Predictions)
    return(mylist)
  }
  else{
  
  Predictions <- predict(Opt.Trees$t.object, XTesting)
  error       <- 1-sum(diag(table(Predictions, as.factor(YTesting))))/sum(table(Predictions, as.factor(YTesting)))
  confusion <- table("Predicted Class"=Predictions, "True Class"=as.factor(YTesting))
  
  mylist <- list("Error.Rate"=error,"Confusion.Matrix"=confusion,"Predicted.Class.Labels"=Predictions, "Trees.Used"=Opt.Trees$t.object$ntree)
  return(mylist)
  }
}
