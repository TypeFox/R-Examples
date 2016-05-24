mWindow_main <-
function (data, dimens, alpha, windowSize) {

  formatted_data <- array (0, c(dim(data)[1], dim(data)[2]))
  for (j in 1:(dim(data)[2]-1)) {
    v <- as.factor(data[,j])
    if (length(levels(v)) > dimens[j]) {
      stop(paste("The dimens vector does not agree with the data. For example, dimens[",j,"] must be increased to at least ", length(levels(v)), ".", sep = ""))
    }
    formatted_data[,j] <- as.double(v) - 1
  }

  v <- as.factor(data[,dim(data)[2]])
    if (length(levels(v)) != 2)
     stop ("Response must be binary")
  formatted_data[,dim(data)[2]] <- as.double(v) - 1
  
  formatted_data <- as.data.frame(formatted_data)
  colnames(formatted_data) <- colnames(data)
  
  data <- formatted_data
  tData <- t(data)
  nVars <- dim(data)[2]
  
  prettyFormulas <- logMargLik <- array (NA, c(nVars - windowSize))
  currentModel <- rep (1:windowSize)  

  for (i in 1:length(prettyFormulas)) {
    prettyFormulas[i] <- prettyFormula (c(currentModel, nVars), colnames(data))
    logMargLik[i] <- findLogMargLik (c(currentModel, nVars), tData, dimens, alpha)   
    currentModel <- currentModel + 1
  }

  masterList <- data.frame(formula = prettyFormulas, logMargLik = logMargLik, stringsAsFactors = F)
  masterList <- masterList[order(masterList$logMargLik, decreasing = T),]
  rownames(masterList) <- rep(1:dim(masterList)[1])
  return(masterList)

}
