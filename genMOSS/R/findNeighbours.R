findNeighbours <-
function (model, potVars, maxVars, confVars) {
  
  model <- na.omit(model)
  nVars <- length(model)  
  modelComp <- which (!(potVars %in% model))
  nVarsComp <- length(modelComp)   
  nConfVars <- length(confVars)   

  neighbourList <- c()

  if (nVars < maxVars) {
    neighbourList <- repRow (model, nVarsComp)
    neighbourList <- cbind(neighbourList, modelComp)
    neighbourList <- cbind(neighbourList, array(NA, c(nVarsComp, maxVars - nVars - 1)))
  }

  for (i in 1:(nVars - nConfVars - 1)) {
      tempModel <- model[-i]
      tempList <- repRow(tempModel, nVarsComp)
      tempList <- cbind(tempList, modelComp, array(NA, c(nVarsComp, maxVars - nVars)))
      neighbourList <- rbind(neighbourList, tempList) 
      if (nVars - nConfVars - 1 > 1) 
        neighbourList <- rbind(neighbourList, c(tempModel, rep(NA, maxVars - nVars + 1)))    
  }
     
  neighbourList <- t(apply(neighbourList, 1, function (x) {sort(x, na.last = T)}))

  return(neighbourList)
 
}
