MOSS_GWAS_main <-
function (alpha, c, cPrime, q, replicates, maxVars, data, dimens, confVars, k) {

  # To prevent the warning from the glm function when fitting a log-linear model for a contigency table 
  # with non-integer values (i.e. the posterior and prior tables) 
  options(warn = -1)
    
  # if (!is.null(k)) {
  #  ROCR <- require(ROCR)
  #  if (!ROCR) 
  #    stop ("Please install the ROCR package")
  # }

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
  nPotVars <- dim(data)[2]
  potVars <- rep(1:nPotVars)
  nConfVars <- length(confVars) 
 
  # Move confounding variables and response to the end of data
  a <- c(which(colnames(data) %in% confVars), nPotVars)
  data <- data.frame(data[,-a, drop = F], data[, a, drop = F])
  dimens <- c(dimens[-a], dimens[a])
  varNames <- colnames(data)
  tData <- t(data)

  masterList <- c()
  modelsInMasterList <- c() 
  cat ("\nsearching for top regression models...\n")

  for (r in 1:replicates) {
  
    # generate a random starting point

    #randomModel <- c(nPotVars, sample(nPotVars-1, sample (maxVars-1, size = 1)))
    randomModel <- c(seq(nPotVars - nConfVars, nPotVars), sample(nPotVars - nConfVars - 1, sample (maxVars - nConfVars - 1, size = 1)))
    nVars <- length(randomModel)
    models <- matrix (nrow = 1, ncol = maxVars)
    models[1,] <- c(randomModel, rep(NA, maxVars - nVars))
    models[1,] <- sort(models[1,], na.last = T)
    prettyFormulas <- prettyFormula (models[1,], varNames)
    logMargLik <- findLogMargLik(models[1,], tData, dimens, alpha)
    explored <- 0
     
    iteration <- 1
    
    cat ("\n")

    while(1) {

      numUnExploredModels <- sum(1 - explored)
      outputMessage1 <- paste ("replicate [", r, "], ", "iteration [", iteration, "].", sep = "")
      outputMessage2 <- paste ("models in list [", length(prettyFormulas), "], ", "not studied [", numUnExploredModels, "].", sep = "")
      cat(outputMessage1, "\n", outputMessage2, "\n\n", sep = "")

      if (numUnExploredModels == 0) {
        currentList <- data.frame(formula = prettyFormulas, logMargLik = logMargLik, stringsAsFactors = F)
        toKeep <- logMargLik >= log(c) + max(logMargLik)
        currentList <- currentList [toKeep,,drop = F]
        models <- models[toKeep,,drop = F]
        masterList <- rbind(masterList, currentList) 
        modelsInMasterList <- rbind(modelsInMasterList, models)
        break 
      }

      unExploredModels <- which(explored == 0)   
 
      # choose a model

      if (length(unExploredModels) == 1) {
        m <- unExploredModels  
      }
      else {
        unExploredmargLik <- logMargLik[explored == 0]
        unExploredmargLik <- unExploredmargLik - max(unExploredmargLik)
        unExploredmargLik <- exp(unExploredmargLik)
        m <- sample (x = unExploredModels, size = 1, prob = unExploredmargLik)
      }

      explored[m] <- 1
      neighbourList <- findNeighbours (models[m,], potVars, maxVars, confVars)   

      # visit all neighbouring models

      toKeep <- vector (mode = "logical", dim(neighbourList)[1])
      currentPrettyFormulas <- vector (mode = "character", dim(neighbourList)[1])

      for (i in 1:dim(neighbourList)[1]) {
        currentPrettyFormulas[i] <- prettyFormula (neighbourList[i,], varNames) 
        toKeep[i] <- !(currentPrettyFormulas[i] %in% prettyFormulas)
      }    
      
      neighbourList <- neighbourList[toKeep,,drop = F]
      
      if (dim(neighbourList)[1] > 0) {
        
        prettyFormulas <- c(prettyFormulas, currentPrettyFormulas[toKeep])  
        currentLogMargLik <- vector (mode = "numeric", dim(neighbourList)[1]) 
        models <- rbind (models, neighbourList)

        for (i in 1:dim(neighbourList)[1]) 
          currentLogMargLik[i] <- findLogMargLik (neighbourList[i,], tData, dimens, alpha) 

        logMargLik <- c(logMargLik, currentLogMargLik)
        explored <- c(explored, array (0,dim(neighbourList)[1]))                     
      
      }

      criteria1 <- logMargLik >= log(cPrime) + max(logMargLik)
      criteria2 <- ifelse(logMargLik >= log(c) + max(logMargLik), 1, rbinom(1,1,1-q)) 
      toKeep <- criteria1 & criteria2
      models <- models[toKeep,,drop = F]        
      prettyFormulas <- prettyFormulas[toKeep]
      logMargLik <- logMargLik[toKeep]   
      explored <- explored[toKeep]
      iteration <- iteration + 1
      
    }

  } 

  masterList <- unique(masterList)
  modelsInMasterList <- unique(modelsInMasterList)

  normMargLik <- masterList$logMargLik - max(masterList$logMargLik)
  normMargLik <- exp(normMargLik)
  normMargLik <- normMargLik/ sum(normMargLik)
      
  vars <- unique(sort(modelsInMasterList))
  #vars <- vars[-length(vars)]
  vars <- vars[1:(length(vars) - nConfVars - 1)]
  nVars <- length(vars)
  postIncProb <- array(0, c(nVars))

  m <- 1
  for (i in vars) {
    toSum <- ceiling(which(t(modelsInMasterList) == i) / maxVars)
    postIncProb[m] <- sum(normMargLik[toSum])
    m <- m + 1
  }

  postIncProbList <- data.frame(variable = varNames[vars], postIncProb = postIncProb, stringsAsFactors = F)

  order <- order(normMargLik, decreasing = T)
  masterList <- masterList[order,,drop = F]
  modelsInMasterList <- modelsInMasterList[order,,drop = F]
  postIncProbList <- postIncProbList[order(postIncProbList$postIncProb, decreasing = T),,drop = F]
  
  rownames(masterList) <- rep(1:dim(masterList)[1])
  rownames(modelsInMasterList) <- rep(1:dim(modelsInMasterList)[1])
  rownames(postIncProbList) <- rep(1:dim(postIncProbList)[1])
  
  if (is.null(k)) {
    x <- list(topRegressions = masterList, postIncProbs = postIncProbList)
  }
  else {
    interactionModels <- c()
  
    cat ("searching for top hierarchical log-linear models...\n")

    fits <- list()

    for (i in 1:dim(modelsInMasterList)[1]) {
      margData <- data[,na.omit(modelsInMasterList[i,]),drop = F]
      for (j in 1:dim(margData)[2])
        margData[,j] <- factor(margData[,j], levels = rep(0:(dimens[modelsInMasterList[i,j]]-1)))
      margData <- as.data.frame.table(table(margData))
      colnames(margData)[dim(margData)[2]] <- "freq"
      interactionModels <- rbind(interactionModels, MOSS_Hierarchical (alpha = alpha, c = 0.1, cPrime = cPrime, q = q, replicates = 5, data = margData)) 
      fits[[i]] <- glm (formula = interactionModels$V2[i], family = poisson(), data = margData) 
    }

    cat("performing cross validation...\n")
    cvDiag <- cvFunc (data, dimens, interactionModels, modelsInMasterList, normMargLik, k, alpha)
    interactionModels$V2 = NULL
    colnames(interactionModels) <- c("formula", "logMargLik")

    x <- list(topRegressions = masterList, postIncProbs = postIncProbList, interactionModels = interactionModels, fits = fits, cvMatrix = cvDiag$cvMatrix, cvDiag = data.frame(acc = cvDiag$acc, tpr = cvDiag$tpr, fpr = cvDiag$fpr , auc = cvDiag$auc))
  
    for (i in 1:dim(x$interactionModels)[1])
      x$fits[[i]]$call <- interactionModels[[i]]  

  }   
  
  options (warn = 0)

  return(x)

}
