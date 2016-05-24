Predict.esknnClass <-
function(optModels, xtest,ytest=NULL,k=NULL)
{ 
  mod <- function(x){
    vt <- table(x)
    as.numeric(names(vt[vt == max(vt)]))}
  
  k <- ifelse(is.null(k),3,k)
  fit<-list()
  zpred<-list()
  zprob<-list()
  ### selecting odd number of models in final ensemble to break th ties if any in voting
  len <- length(optModels$fsfinal)
  if(len%%2==0)
    len <- len-1  
  ### predicting test data using final ensemble
  for(z in 1:len) 
  {
    fit<- knn3Train(optModels$trainfinal[[z]][,names(optModels$trainfinal[[z]])!="Class"],xtest<-xtest[,optModels$fsfinal[[z]]],optModels$trainfinal[[z]]$Class,k=k)
    ## extract class vector and probability vector from knn3
    zpred[[z]]<-as.factor(fit[1:length(fit)])                    ### class labels
  }
  ##   binding selected z models 
  mclass<- do.call("cbind",zpred)
  predClass<-apply(mclass,1,mod)
  ##
  ## Calculating classification Error
  ##
  if(is.null(ytest)){
    return(list("PredClass"= predClass))
  }
  else{
  conf=table("True.Class"=as.numeric(ytest)-1,"Predicted.Class"=as.numeric(predClass)-1)   
  err=1-(sum(diag(conf))/nrow(xtest))
  return(list("PredClass"= predClass,"ConfMatrix"=conf,"ClassError"=err))
}
}
