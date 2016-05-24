cvFunc <-
function (data, dimens, interactionModels, modelsInMasterList, regPostProb, k, alpha) {
  
  newOrder <- sample (x = rep(1:dim(data)[1]), size = dim(data)[1], replace = F)
  data <- data[newOrder,]   
  partitionSize <- dim(data)[1] %/% k
  confusionMatrix <- list()

  for (i in 1:k) {

    first <- 1 + (i-1) * partitionSize
    if (i == k)
      last <- dim(data)[1]
    else
      last <- i * partitionSize
    test.indices <- rep(first:last)
    test.data <- data[test.indices,,drop = F]
    train.data <- data[-test.indices,,drop = F] 
    fitted.values <- fitInteractionModels (train.data, dimens, interactionModels, modelsInMasterList, alpha)
    diseaseProb <- array (0, dim(test.data)[1])
    prediction <- array (0, dim(test.data)[1])
    for (r in 1:dim(test.data)[1]) {
      for (s in 1:dim(modelsInMasterList)[1]) {
        obsVector <- test.data[r, na.omit(modelsInMasterList[s,])]
        obsVector[length(obsVector)] <- 0
        v0 <- fitted.values[[s]][findIndex(obsVector, dimens[na.omit(modelsInMasterList[s,])])]
        obsVector[length(obsVector)] <- 1
        v1 <- fitted.values[[s]][findIndex(obsVector, dimens[na.omit(modelsInMasterList[s,])])]
        diseaseProb[r] <- diseaseProb[r] + regPostProb[s] * v1 / (v0 + v1)
      }
      prediction[r] <- round(diseaseProb[r])
    }
    actual <- test.data[,dim(test.data)[2]]
    actual <- factor(actual, levels = c(0,1))
    prediction <- factor(prediction, levels = c(0,1))
    confusionMatrix[[i]] <- data.frame(prediciton = prediction, actual = actual)
    confusionMatrix[[i]] <- as.data.frame.table(table(confusionMatrix[[i]]))
    colnames(confusionMatrix[[i]])[dim(confusionMatrix[[i]])[2]] <- "prop"
    confusionMatrix[[i]]$prop <- confusionMatrix[[i]]$prop / dim(test.data)[1]
  }

  avgConfusionMatrix <- confusionMatrix[[1]]

  for (i in 2:k) 
    avgConfusionMatrix$prop <- avgConfusionMatrix$prop + confusionMatrix[[i]]$prop
         
  avgConfusionMatrix$prop <- avgConfusionMatrix$prop / k

  return(avgConfusionMatrix)

}
