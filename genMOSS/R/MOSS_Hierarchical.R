MOSS_Hierarchical <-
function (startList = NULL, p = 0.2, alpha = 1, c = 0.1, cPrime = 0.0001, q = 0.1, replicates = 5, data) {
  
  tools <- list()

  varNames <- colnames(data)[which(colnames(data)!= "freq")]
  n <- length(varNames)
  varSets <- decToBin (0:(2**n-1),n)
  colnames(varSets) <- varNames
  lenVarSets <- rowSums(varSets)
  nVarSets <- 2 ** n   

  # lattice
  
  downLinks <- array(NA, c(nVarSets,n))
  nDownLinks <- lenVarSets
  upLinks <- array(NA,c(nVarSets,n))
  nUpLinks <- n - lenVarSets

  # downLinks

  for(i in 1:nVarSets) {
    k = 1
    for(j in 1:n) {
      if(varSets[i,j] == 1) {
        varSets[i,j] <- 0
        downLinks[i,k] <- binToDec(varSets[i,]) + 1
        k <- k + 1
        varSets[i,j] <- 1
      }
    }
  }

  # upLinks

  for(i in 1:nVarSets) {
    k = 1
    for(j in 1:n) {
      if(varSets[i,j] == 0) {
        varSets[i,j] <- 1
        upLinks[i,k] <- binToDec(varSets[i,]) + 1
        k <- k + 1
        varSets[i,j] <- 0
      }
    }
  }  

  tools <- list(varNames = varNames, n = n, varSets = varSets, lenVarSets = lenVarSets, nVarSets = nVarSets, downLinks = downLinks, nDownLinks = nDownLinks, upLinks = upLinks, nUpLinks = nUpLinks)

  postData <- priorData <- data
  postData$freq <- data$freq + alpha / length(data$freq)
  priorData$freq <- array (alpha / length(data$freq), c(length(data$freq))) 

  sizeOfStartList <- length(startList)
  
  masterList <- c()

  for (r in 1:replicates) {
  
    models <- matrix (nrow = 1, ncol = nVarSets)
    generators <- matrix (nrow = 1, ncol = nVarSets)
    dualGenerators <- matrix (nrow = 1, ncol = nVarSets)
    models[1,] <- randomHierModel(p, tools)
    generators[1,] <- findGenerators(models[1,], tools)
    dualGenerators[1,] <- findDualGenerators(models[1,], tools)
    formulas <- findFormula(generators[1,], tools)
    logMargLik <- logLaplace(formulas, postData, tools) - logLaplace(formulas, priorData, tools) 
    explored <- 0    
          
    iteration <- 1
    #cat ("\n")

    while(1) {

      numUnExploredModels <- sum(1 - explored)
      #outputMessage1 <- paste ("replicate [", r, "], ", "iteration [", iteration, "].", sep = "")
      #outputMessage2 <- paste ("models in list [", length(formulas), "], ", "not studied [", numUnExploredModels, "].", sep = "")
      #cat(outputMessage1, "\n", outputMessage2, "\n\n", sep = "")

      if (sum(explored) == length(explored)) {
        prettyHierFormulas <- vector(mode = "character", length(formulas))
        for (i in 1:length(prettyHierFormulas))
          prettyHierFormulas[i] <- prettyHierFormula (generators[i,], tools)
        currentList <- data.frame(V1 = prettyHierFormulas, V2 = formulas, V3 = logMargLik, stringsAsFactors = F)
        currentList <- currentList [logMargLik >= log(c) + max(logMargLik),]
        masterList <- rbind(masterList, currentList) 
        break 
      }

      unExploredModels <- which(explored == 0)   
 
      if (length(unExploredModels) == 1) {
        m <- unExploredModels  
      }
      else {
        unExploredPostProb <- logMargLik[explored == 0]
        unExploredPostProb <- unExploredPostProb - max(unExploredPostProb)
        unExploredPostProb <- exp(unExploredPostProb)
        m <- sample (x = unExploredModels, size = 1, prob = unExploredPostProb)
      }

      explored[m] <- 1
      neighbourList <- findHierNeighbours (models[m,], generators[m,], dualGenerators[m,], tools)  
    
      for (i in 1:dim(neighbourList)[1]) {
 
        currentNeighbourGenerators <- findGenerators (neighbourList[i,], tools)
        currentNeighbourFormula <- findFormula(currentNeighbourGenerators, tools)
        inList <- currentNeighbourFormula %in% formulas
                
        if (!inList) {
          models <- rbind(models, neighbourList[i,])
          formulas <- c(formulas, currentNeighbourFormula)
          generators <- rbind(generators, currentNeighbourGenerators)
          currentNeighbourDualGenerators <- findDualGenerators(currentNeighbourGenerators, tools)
          dualGenerators <- rbind(dualGenerators, currentNeighbourDualGenerators)
          logMargLik <- c(logMargLik, logLaplace(as.formula(currentNeighbourFormula), postData, tools) - logLaplace(as.formula(currentNeighbourFormula), priorData, tools))
          explored <- c(explored, 0) 
        }
      }
 
      criteria1 <- logMargLik >= log(cPrime) + max(logMargLik)
      criteria2 <- ifelse(logMargLik >= log(c) + max(logMargLik), 1, rbinom(1,1,1-q)) 
      toKeep <- criteria1 & criteria2
      models <- models[toKeep,,drop = F]        
      formulas <- formulas[toKeep]
      generators <- generators[toKeep,,drop = F]
      dualGenerators <- dualGenerators[toKeep,,drop = F]
      logMargLik <- logMargLik[toKeep]   
      explored <- explored[toKeep]
      iteration <- iteration + 1
      
    }
  
    explored <- rep(0,length(explored))
    
  }  

  masterList <- unique(masterList)
  masterList <- masterList[order(masterList$V3, decreasing = T),]
  row.names(masterList) <- rep(1:dim(masterList)[1])

  return(masterList[1,,drop = F])  
 
}
