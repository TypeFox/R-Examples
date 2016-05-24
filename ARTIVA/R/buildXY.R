buildXY <-
function(targetData, parentData, dataDescription, m, dyn){
  ### Build response Y and predictor X
  if(is.null(dataDescription)){
    Y=targetData
    X=parentData
  }else{
    pos= order(dataDescription)
    Y=targetData[pos]
    X=parentData[, pos]
  }

  lgth=length(targetData)
  Y=Y[(m*dyn+1):lgth]
  ##Y=Y[-c(1:(m*dyn))]

 
  ## updated by Sophie 07/07/2011
  if(dim(X)[1]==1){
    ## updated by Sophie 12/07/2011
    X=matrix(X[,1:(lgth-(m*dyn))],1,length(Y))
  }else{
    ## updated by Sophie 12/07/2011
    X=X[,1:(lgth-(m*dyn))]
  }

  X=t(X)

  ## add a constant vector to predictor data
 
  X = cbind(X,array(1,dim(X)[1]))
  # return formatted data
  return(list(X=X,Y=Y))
}
