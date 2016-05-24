Predict.OTProb <-
function(Opt.Trees,XTesting, YTesting=NULL){
  
  Pr=predict(Opt.Trees$t.object, XTesting,type='prob')[,2]
  
  if(is.null(YTesting)){
    mylist <- list("Estimated.Probabilities"=Pr)
    return(mylist)
  }
  else{
  BS=mean((as.vector(Pr)-(as.numeric(YTesting)))^2)
  
  mylist <- list("Brier.Score"=BS,"Estimated.Probabilities"=Pr)
  return(mylist)
  }
}
