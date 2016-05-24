runBackShift <- function(X, environment, parentsOf, variableSelMat, pointEst, 
                         setOptions, verbose, result){
  
  # additional options for backShift
  optionsList <- list("covariance"=TRUE, "ev"=0, "threshold"=0.75, "nsim"=100, 
                      "sampleSettings"=1/sqrt(2), "sampleObservations"=1/sqrt(2),
                      "nodewise"=TRUE, "tolerance"=10^(-4), "baseSettingEnv"=1)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  p <- ncol(X)
  
  if(nrow(X) < p) 
    stop("backShift not suitable if there are more variables 
         than observations")
  if(!is.null(variableSelMat)) 
    warning("option 'variableSelMat' not implemented for 
            'backShift' -- using all variables")
  
  res <- try(backShift::backShift(
     X, environment, covariance=optionsList$covariance, 
     ev=optionsList$ev, threshold=optionsList$threshold, 
     nsim=optionsList$nsim, 
     sampleSettings=optionsList$sampleSettings, 
     sampleObservations=optionsList$sampleObservations, 
     nodewise=optionsList$nodewise, 
     tolerance=optionsList$tolerance,
     baseSettingEnv = optionsList$baseSettingEnv, 
     verbose = verbose), 
     silent = FALSE)
  
  if(inherits(res, "try-error")){
    warning("backShift -- no stable model found. Possible model 
            mispecification. Returning the empty graph.\n")
    
    res<- list(Ahat=0*diag(p), AhatAdjacency = 0*diag(p), 
               varianceEnv = matrix(0, nrow = length(unique(environment)), 
                                    ncol = p))
  }
  for (k in 1:length(parentsOf)){
    if(optionsList$ev == 0)
      result[[k]] <- (wh <- which(res$Ahat[, k]!=0))
    else
      result[[k]] <- (wh <- which(res$AhatAdjacency[, k]!=0))
    
    if(pointEst) attr(result[[k]],"coefficients") <- res$Ahat[ wh,k ]
  }
  
  result
}