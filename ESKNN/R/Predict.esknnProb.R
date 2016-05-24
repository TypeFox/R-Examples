Predict.esknnProb <-
function(optModels, xtest,ytest=NULL,k=NULL)
{ 
  k <- ifelse(is.null(k),3,k)  
  fit<-list()
  zprob<-list()
  ### selecting odd number of models in final ensemble to break th ties if any in voting
  len <- length(optModels$fsfinal)
  if(len%%2==0)
    len <- len-1
  
  ### predicting test data using final ensemble

  for(z in 1:len)
    
  {
    fit<- knn3Train(optModels$trainfinal[[z]][,names(optModels$trainfinal[[z]])!="Class"],xtest<-xtest[,optModels$fsfinal[[z]]],optModels$trainfinal[[z]]$Class,k=k)
    ## extract  probability vector from knn3
    zprob[[z]]<-attributes(fit)$prob[,2]                         ### class probabilities
  }
  ##   binding selected z models 
  mprob<-  do.call("cbind",zprob)
  predProb<- apply(mprob,1,mean)
  predProb<-as.vector(predProb)  
  ##
  ###Computing Brier Score
  ##
  
  if(is.null(ytest)){
    return(list("PredProb"=predProb))
  }
  else{
  ytest<-as.numeric(as.factor(ytest))-1
  BS=mean((predProb-ytest)^2)
  ##
  ##Returning List
  ##
  return(list("PredProb"=predProb, "BrierScore"=BS))
}
}
