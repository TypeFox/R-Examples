runRegression <- function(X, parentsOf, variableSelMat, pointEst, setOptions, 
                          verbose, result){
  
  # additional options for regression
  optionsList <- list("selfselect"=NULL)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  for (k in 1:length(parentsOf)){
    if(round(k/100)==(k/100)) cat(" ",k)
    
    possibleVar <- (1:ncol(X))
    
    possibleVar <- (1:ncol(X))
    removeVar <- parentsOf[k]
    if(!is.null(variableSelMat) & is.null(optionsList$selfselect)){
      if(ncol(variableSelMat)==length(parentsOf)) 
        selc <- k 
      else 
        selc <- which( (1:ncol(X)) == parentsOf[k])
      removeVar <- 
        unique(c(removeVar,which( !variableSelMat[,selc] )))
    }else{
      if(!is.null(optionsList$selfselect)){
        gl <- glmnet::glmnet(X[, possibleVar[-removeVar],drop=FALSE], 
                             as.numeric(X[, parentsOf[k]]))
        nnz <- apply(coef(gl)!=0, 2,sum)
        beta <- 
          coef(gl, s=gl$lambda[sum(nnz<=optionsList$selfselect)])[-1]
        removeVar <- 
          c(removeVar,(possibleVar[-removeVar])[which(beta==0)])
      }
    }
    if(length(removeVar)>0) 
      possibleVar <- possibleVar[-removeVar]
    
    
    parents <- numeric(0)
    if(length(possibleVar)>1){
      gl <- glmnet::cv.glmnet(X[ ,possibleVar,drop=FALSE], 
                              as.numeric( X[,parentsOf[k]]),
                              intercept=TRUE)
      beta <- as.numeric(coef(gl))[-1]
      parents <- possibleVar[which(beta!=0)]
    }else{
      beta <- 0
    }
    result[[k]] <- parents
    if(pointEst) 
      attr(result[[k]],"coefficients") <- beta[beta!=0]
    
  }
  
  result
}