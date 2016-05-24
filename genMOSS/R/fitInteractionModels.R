fitInteractionModels <-
function (train.data, dimens, interactionModels, modelsInMasterList, alpha) {

  fitted.values <- list()

  for (i in 1:dim(interactionModels)[1]) {
    
    margData <- train.data[,na.omit(modelsInMasterList[i,]),drop = F]
    for (j in 1:dim(margData)[2])
      margData[,j] <- factor(margData[,j], levels = rep(0:(dimens[modelsInMasterList[i,j]]-1)))
    margData <- as.data.frame.table(table(margData))
    colnames(margData)[dim(margData)[2]] <- "freq"

    margData$freq <- margData$freq + alpha / length(margData$freq) 
    fitted.values[[i]] <- glm (formula = interactionModels$V2[i], family = poisson(), data = margData)$fitted.values
    
  }
  
  return(fitted.values)

}
