Predict.OTReg <-
function(Opt.Trees,XTesting, YTesting){
  
  
  Predictions=predict(Opt.Trees$t.object, XTesting)
  if(is.null(YTesting)){
    mylist <- list("Pr.Values"=Predictions)
    return(mylist)
  }
  
  else{
  SST <- sum((YTesting-mean(YTesting))^2)
  SSE <- sum((YTesting-Predictions)^2)
  R.Squared=1-SSE/SST
  Unexplained <- SSE/SST #R-squared=1-SSE/SST
  
  mylist <- list("R.Squared"=R.Squared,"Unexp.Variations"=Unexplained,"Pr.Values"=Predictions, "Trees.Used"=Opt.Trees$t.object$ntree)
  return(mylist)
  }
}
