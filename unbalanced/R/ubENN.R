ubENN <-
function(X,Y, k=3, verbose=TRUE){
  
  stopifnot(k > 0, class(verbose) == "logical", all(unique(Y) %in% c(0, 1)))
  
  #only numeric features are allowed
  if(any(sapply(X,is.numeric)==FALSE))
    stop("only numeric features are allowed to compute nearest neighbors")
  
  i.1<-which(Y==1)
  i.0<-which(Y==0)
  
  if(length(i.0)==0){
    #if there are no 0 obs then don't do anything
    if(verbose) 
      cat("Warning: No negative instances \n")
    return(list(X=X,Y=Y))
  }
  
  #removes only example from the majority class
  timeRemove<-system.time({
    out.hat <- FNN::knn(train=X,test=X[i.0,], cl=Y, k=k+1,prob=TRUE) #the 1-nn is the point itself therefore we need k+1
    proba.hat <- attr(out.hat, "prob")
    levels(out.hat) <- c(0,1)
    prob.th <- k/(k+1)
    id.miss <- which((Y[i.0]!=out.hat) & (proba.hat>=prob.th))
  })	 
  if(verbose) 
    cat("Number of instances removed from majority class with ENN:",length(id.miss),
        "\t Time needed:",round(as.numeric(timeRemove["elapsed"]),digits=2),"\n")
  
  if(length(id.miss)==0) 
    id.keep.0<-i.0
  else 
    id.keep.0<-setdiff(i.0,i.0[id.miss])
  
  Id<-c(id.keep.0,i.1)
  Id<-sort(Id)
  id.removed<-setdiff(1:nrow(X),Id)
  
  if (is.vector(X)!=TRUE) X=X[Id,]
  else  X=X[Id]
  Y=Y[Id]
  
  return(list(X=X,Y=Y,id.rm=id.removed))
}
