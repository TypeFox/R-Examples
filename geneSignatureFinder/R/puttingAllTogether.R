puttingAllTogether <-
function (workingFile = "", nchips = 2, alpha = 0.05, saveSurvivalCurvesPlot = FALSE, 
            saveIndividualSignature = FALSE, saveImportancePlot = FALSE) 
  {
    if (areDataNotLoaded()) 
      return(NULL)
    
    if (workingFile != "") 
      message(paste("Using", workingFile, "data."))
    n <- nrow(geData)
    m <- ncol(geData)
    message(paste("Found", n, "samples and", m, "gene expression levels."))
    genesToUse <- !is.na(geData)
    genesToUse <- apply(genesToUse, 2, sum)/n
    toRemove <- which(genesToUse <= 0.75)
    genesToUse <- which(genesToUse > 0.75)
    if (length(genesToUse) != m) {
      message(paste(m - length(genesToUse), "genes are going to be removed from the data matrix (more than 75% of missing values)."))
      message(paste("Removing:", paste(names(toRemove), collapse = ", ")))
      geData <- geData[, genesToUse]
      m <- ncol(geData)
    }
    message(paste("Working on", n, "samples and", m, "gene expression levels."))
    message("Starting the scanning of the genes (seeds finder step)")
    aMakeCluster <- makeCluster(nchips)
    aSeedsFinder <- seedsFinder(cpuCluster = aMakeCluster, cutoff = 1.95)
    attr(aSeedsFinder, "Creation date") <- date()
    
    
    fName <- paste(workingFile, "AnsSeedsFinder.RData", sep = "")
    message("Saving the result of the scanning in ", fName, " of the working directory.")
    save(aSeedsFinder, workingFile, file = fName)
    #####################################
    message("sum(tmp <- is.na(aSeedsFinder[, 'bic1'])) ", sum(tmp <- is.na(aSeedsFinder[,"bic1"])))
    if(sum(tmp) > 0) aSeedsFinder[which(tmp), "bic1"] <- -Inf
    
    message("sum(tmp <- is.na(aSeedsFinder[, 'bic2'])) ", sum(tmp <- is.na(aSeedsFinder[,"bic2"])))
    if(sum(tmp) > 0) aSeedsFinder[which(tmp),"bic2"] <- Inf
    
    
    bimodalGenes <- aSeedsFinder[, "bic1"] > aSeedsFinder[, "bic2"]
    if (sum(bimodalGenes, na.rm = TRUE) < 1) {
      message("No gene have been found bi-modal")
      bimodalGenes <- rep(TRUE, n)
    }
    else message(paste("The ", round(sum(bimodalGenes, na.rm = TRUE)/m, 
                                     4) * 100, "% of the genes have been found bi-modal", 
                       sep = ""))
    
    significantGenes <- BHcorrection(aSeedsFinder[, "pValue"]) < alpha
    if (sum(significantGenes) == 0) {
      message("The B & H correction of the p-values returns no significant genes.")
      significantGenes <- aSeedsFinder[, "pValue"] < alpha
      if (sum(significantGenes) == 0) {
        message("No significant genes without p-value correction at level.", 
                round(alpha, 4), "\n")
        message("The minimum of the p-values is ", round(min(aSeedsFinder[, 
                                                                          "pValue"]), 4), "\n")
        message("The procedure stops.")
        return(NULL)
      }
      message(paste("The analysis proceds on the ", round(sum(significantGenes)/m, 
                                                          4) * 100, "% of the genes found significant (without correction)", 
                    sep = ""))
    }
    else message(paste("The analysis proceds on the ", round(sum(significantGenes)/m, 
                                                             4) * 100, "% of the genes found significant (with correction)", 
                       sep = ""))
    if (sum(bimodalGenes, na.rm = TRUE) == n) {
      seedGenes <- significantGenes
      message("The analysis proceds on the significant genes (not bi-modal).")
    } else {
      seedGenes <- significantGenes * bimodalGenes
      message(paste("The ", round(sum(seedGenes, na.rm = TRUE)/m, 
                                  4) * 100, "% of the genes  has been found significant and bimodal.", 
                    sep = ""))
    }
    seedGenes <- which(seedGenes == 1)
    seedGenesNames <- names(seedGenes)
    message(paste("The seed-genes are:", paste(seedGenesNames, 
                                               collapse = ", ")))
    message("Starting the development of the signature from all seed-genes.")
    K <- length(seedGenes)
    searchResults <- vector("list", K)
    names(searchResults) <- seedGenesNames
    
    aMakeCluster <- makeCluster(nchips)
    for (k in 1:K) {
      message(paste("Developing from the seed:", seedGenesNames[k]))
      if (workingFile != "") 
        fName <- paste(seedGenesNames[k], "@", workingFile, 
                       sep = "")
      else fName <- paste(seedGenesNames[k], sep = "")
      
      aSignatureFinder <- signatureFinder(seedGenes[k], cpuCluster = aMakeCluster, stopCpuCluster = FALSE)
      
      if (saveSurvivalCurvesPlot) {
        sf <- survfit(stData ~ aSignatureFinder$classification)
        plot(aSignatureFinder)
        dev.copy2pdf(file = paste(fName, "SurvivalCurves.pdf", sep = ""))
      }
      
      if (length(aSignatureFinder$signature) > 1) {
        aSignatureFinder <- importance(aSignatureFinder, cpuCluster = aMakeCluster, stopCpuCluster = FALSE)
        
        if (saveImportancePlot) {
          barplot(aSignatureFinder$importance, main = "Importance based on L1GeneOut", 
                  sub = paste("Signature starting from:", aSignatureFinder$startingSignature))
          dev.copy2pdf(file = paste(fName, "Importance.pdf", sep = ""))
        }
      } else message("Importance and test are not computed on this signature because of length = 1.")
      
      aSignatureFinder$workingFile <- workingFile
      if (saveIndividualSignature) 
        save(aSignatureFinder, file = paste(fName, ".RData", sep = ""))
      
      searchResults[[k]] <- aSignatureFinder
    }
    
    stopCluster(aMakeCluster)
    
    fName <- paste(workingFile, "SearchResults.RData", sep = "")
    message("Saving the the signatures found in ", fName, " of the working directory.")
    
    save(searchResults, seedGenesNames, K, workingFile, file = fName)
    ###############################################
    signaturesTable <- searchResultsSummaryTable(searchResults)
    ensemble <- ensembleTable(searchResults)
    #############################  
    write.csv(signaturesTable, file = paste(workingFile, "SignatureTable.csv", sep = ""))
    write.csv(ensemble, file = paste(workingFile, "Ensemble.csv", sep = ""))
    message("The procedure stops.")
    return(NULL)
}
