runLINGAM <- function(X, parentsOf, pointEst, setOptions, directed, verbose, 
                    result){
  
  # additional options for LINGAM
  # TODO: None?
  
  # adjust according to setOptions if necessary
  # optionsList <- adjustOptions(availableOptions = optionsList, 
                               # optionsToSet = setOptions)
                            
  # variableSelMat not implemented                          
  
  if(nrow(X)<=ncol(X)) 
    stop("LINGAM not suitable for high-dimensional data; 
         need nrow(X) > ncol(X)")
  res <- pcalg::LINGAM(X, verbose=verbose)
  lingammat <- res$Adj
  if(directed) lingammat <- lingammat * (t(lingammat)==0)
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- (wh <- which(lingammat[, parentsOf[k]] == 1)) 
    if(pointEst) 
      attr(result[[k]],"coefficients") <- t(res$B)[ wh,parentsOf[k]]
  }
  
  result
}