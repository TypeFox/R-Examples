runBivariateCAM<- function(X, parentsOf, variableSelMat, pointEst, verbose, 
                           result){
  
  bivcammat <- bivariateCAM(X, parentsOf=parentsOf, 
                            variableSelMat=variableSelMat, 
                            silent = !verbose)
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- (wh <- which(bivcammat$causalParents[, parentsOf[k]]>0))
    if(pointEst)
      attr(result[[k]],"coefficients") <- bivcammat$scoreMat[ wh,parentsOf[k] ]
  }
  
  result
}