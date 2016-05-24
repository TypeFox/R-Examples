cvFunc <-
function (data, dimens, interactionModels, modelsInMasterList, weights, k, alpha) {
  
  newOrder <- sample (x = rep(1:dim(data)[1]), size = dim(data)[1], replace = F)
  data <- data[newOrder,]   
  partitionSize <- dim(data)[1] %/% k
  diseaseProbs <- decisions <- rep (0, c(dim(data)[1]))  

  for (i in 1:k) {

    first <- 1 + (i-1) * partitionSize
    if (i == k)
      last <- dim(data)[1]
    else
      last <- i * partitionSize
    test.indices <- rep(first:last)
    train.data <- data[-test.indices,,drop = F] 
    fitted.values <- fitInteractionModels (train.data, dimens, interactionModels, modelsInMasterList, alpha)
    for (r in test.indices) {
      for (s in 1:dim(modelsInMasterList)[1]) {
        obsVector <- data[r, na.omit(modelsInMasterList[s,])]
        obsVector[length(obsVector)] <- 0
        v0 <- fitted.values[[s]][findIndex(obsVector, dimens[na.omit(modelsInMasterList[s,])])]
        obsVector[length(obsVector)] <- 1
        v1 <- fitted.values[[s]][findIndex(obsVector, dimens[na.omit(modelsInMasterList[s,])])]
        diseaseProbs[r] <- diseaseProbs[r] + weights[s] * v1 / (v0 + v1)
      }
      decisions[r] <- round(diseaseProbs[r])
    }
  }

  phenos <- data[,dim(data)[2]]
  phenos <- factor(phenos, levels = c(0,1))
  decisions <- factor(decisions, levels = c(0,1))
  cvMatrix <- table(pheno = phenos, decision = decisions)

  acc <- (cvMatrix[1] + cvMatrix[4]) / sum(cvMatrix) * 100 
  tpr <- cvMatrix[4] / (cvMatrix[2] + cvMatrix[4]) * 100
  fpr <- cvMatrix[3] / (cvMatrix[1] + cvMatrix[3]) * 100
  auc <- performance (prediction (diseaseProbs, phenos), "auc")
  auc <- as.numeric(attr (auc, "y.values")) * 100

  return (list (cvMatrix = cvMatrix, acc = round(acc, digits = 1), tpr = round(tpr, digits = 1), fpr = round(fpr, digits = 1), auc = round(auc, digits = 1)))
 
}
