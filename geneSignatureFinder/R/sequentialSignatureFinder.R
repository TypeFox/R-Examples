sequentialSignatureFinder <-
function(startingGene, logFileName = "", coeffMissingAllowed = 0.75,
                            subsetToUse = 1:ncol(geData)) {
  n <- nrow(geData)
  m <- ncol(geData)

  toExplore <- rep(FALSE, m)
  toExplore[subsetToUse]  <- TRUE
  toExplore[startingGene]  <- FALSE
  
  logFileName <- paste(logFileName, "SignatureFinderLog.txt", sep = "")
  cat(paste("Working on ", m, " genes observed on ", n, " samples.\n", sep = ""), file = logFileName, append = FALSE)
  cat(paste("Starting on ", ttime <- Sys.time(), "\n", sep = ""), file = logFileName, append = TRUE)
  cat(paste("Starting signature: ", paste(colnames(geData)[startingGene], collapse = ", "), ";\n", sep = ""),
      file = logFileName, append = TRUE)
  
  aClassify <- rep(NA, n)
  ssignature <- startingGene
  
  runs <- 0 
  result <- list()  
  result$signatureName <-  (colnames(geData)[startingGene])[1]
  result$startingSignature <-  colnames(geData)[startingGene]
  result$coeffMissingAllowed <- coeffMissingAllowed

######## da rivedere; si attiva nel caso lo startingGene sia una sequenza
  if(length(ssignature) > 1) 
    notMissing <- apply(!is.na(geData[, ssignature]), 1, sum) else
  notMissing <- !is.na(geData[, ssignature]) + 0
  notMissing <- notMissing > 0

  clusters <- classify(geData[notMissing, ssignature])$clusters   
  tmp1 <- min(survfit(stData[clusters == 1]~ 1)$surv)
  tmp2 <- min(survfit(stData[clusters == 2]~ 1)$surv)
  if(tmp1 > tmp2) {  
    clusters[clusters == 1] <- 0
    clusters[clusters == 2] <- 1
  } else clusters[clusters == 2] <- 0

  runningDistance <- survdiff(stData[notMissing] ~ clusters)$chisq # variato

  result$startingTValue <- runningDistance
  result$startingPValue <- 1 - pchisq(runningDistance, df = 1)
  tmpClassification <- rep(NA, n)
  tmpClassification[notMissing] <- clusters 
  result$startingClassification <- as.factor(tmpClassification)
  levels(result$startingClassification) <- c("good", "poor")

  cat(paste("tValue = ", runningDistance, " (", result$startingPValue, ")\n", sep = ""), file = logFileName, append = TRUE)
####################### end: da rivedere

  #### MAIN LOOP #######################
  exitFromMain <- FALSE
  repeat {
    if(exitFromMain) break
    runs <- runs + 1
    
    
    
    distances <- rep(0, m) #variato
    j <- 1
    while(j <= m) {
      if(!toExplore[j]) { j <- j + 1; next }
      distances[j] <- tValueFun(c(ssignature, j), coeffMissingAllowed = coeffMissingAllowed)
      j <- j + 1
    }
    
    ### INNER LOOP ##########################
    exitFromInner <-  FALSE
    repeat {
      if(exitFromInner) break
      maxDistance <- max(distances) #variato
      if(maxDistance >= runningDistance) { #variato
        candidate <- which(distances == maxDistance)


        if(length(candidate) == 1) {
                                        # qua si controlla se l'aggiunta dei candidati porta a gruppi degeneri
          notMissing <- apply(!is.na(geData[, c(ssignature, candidate)]), 1, sum)
          notMissing <- notMissing > floor(length(c(ssignature, candidate))^coeffMissingAllowed)
          clusters <- classify(geData[notMissing, c(ssignature, candidate)])$clusters
          checkOne <- (min(table(clusters)) > floor(0.1 * n))
          # qua si controlla che le curve di sopravvivenza non si incrocino
          sf0 <- survfit(stData[clusters == 1] ~ 1)$surv
          sf1 <- survfit(stData[clusters == 2] ~ 1)$surv
          checkTwo <- sum(fivenum(sf0) > fivenum(sf1))
          checkTwo <- ((checkTwo == 0) | (checkTwo == 5))
          
          if(checkOne & checkTwo) {
            ssignature <- c(ssignature, candidate)
            toExplore[ssignature] <- FALSE
            runningDistance <- maxDistance
            cat(paste("... improved signature: ",
                      paste(colnames(geData)[ssignature], collapse = ", "),
                      ";\ntValue = ", runningDistance,  " (",
                      1 - pchisq(runningDistance, df = 1), ")\n", sep = ""),
                file = logFileName, append = TRUE)
            exitFromInner <- TRUE
          }
        } else {
          if(length(candidate) > 0.01*m) {
            exitFromInner <- TRUE
            exitFromMain <- TRUE
            break
          }
          tmpCandidate <- findBest(runningDistance, ssignature, candidate)
	  if(length(tmpCandidate) > 0) { 
	    tmp <- survdiff(stData ~ classify(geData[, c(ssignature, tmpCandidate)])$clusters)$chisq
	    if(tmp > runningDistance) {
                                        #qua si dovrebbero aggiungere un po' di controlli
              candidate <- tmpCandidate
              ssignature <- c(ssignature, candidate)
              runningDistance <- tmp
              toExplore[ssignature] <- FALSE
              cat(paste("... improved signature: ",
                        paste(colnames(geData)[ssignature], collapse = ", "),
                        ";\ntValue = ", runningDistance,  " (",
                        1 - pchisq(runningDistance, df = 1), ")\n", sep = ""),
                  file = logFileName, append = TRUE)
              exitFromInner <- TRUE
            }  
          }       
        }
        distances[candidate] <- 0 #variato        
      } else {
        exitFromInner <- TRUE
        exitFromMain <- TRUE
      }
    } ### END of INNER LOOP ###########
  } ### END of MAIN LOOP ###########


  if(length(ssignature) > 1) {
    notMissing <- apply(!is.na(geData[, ssignature]), 1, sum)
    notMissing <- notMissing > floor(length(ssignature)^coeffMissingAllowed)
  } else {
    notMissing <- !is.na(geData[, ssignature]) + 0
    notMissing <- notMissing > 0
  }

  clusters <- rep(NA, n)
  clusters[notMissing] <- classify(geData[notMissing, ssignature])$clusters
  
  clusters <- goodAndPoorClassification(clusters) # 08/04/2012
  #    tmp1 <- min(survfit(stData[clusters[notMissing] == 1]~ 1)$surv)
  #    tmp2 <- min(survfit(stData[clusters[notMissing] == 2]~ 1)$surv)
  #    if(tmp1 > tmp2) {  
  #      clusters[notMissing][clusters[notMissing] == 1] <- 0
  #      clusters[notMissing][clusters[notMissing] == 2] <- 1
  #    } else clusters[notMissing][clusters[notMissing] == 2] <- 0
  
  K <- length(ssignature)
  result$signature <-  colnames(geData)[ssignature]
  result$tValue <- survdiff(stData[notMissing] ~ clusters[notMissing])$chisq 
  result$pValue <- 1-pchisq(result$tValue, df = 1)
  result$signatureIDs <- ssignature
  names(result$signatureIDs) <- result$signature
  result$classification <- clusters#08/04/2012
  #result$classification <- as.factor(clusters)
  #levels(result$classification) <- c("good", "poor")
  cat(paste("\n\nfinal signature: ", paste(colnames(geData)[ssignature], collapse = " "), sep = ""), file = logFileName, append = TRUE)
  cat(paste("\ntValue = ", runningDistance, " (", 1 - pchisq(runningDistance, df = 1), ")\n",  sep = ""), file = logFileName, append = TRUE)
  cat(paste("\nlength of the signature = ", length(ssignature), sep = ""), file = logFileName, append = TRUE)
  cat(paste("\nnumber of joint missing values = ", sum(!notMissing), " (", 100*round(sum(!notMissing)/n,2), "%)", sep = ""), file = logFileName, append = TRUE)
  cat(paste("\n\nEnd of computation at ", t2 <- Sys.time(), "; elapsed time: ", t2 - ttime, ".", sep = ""), file = logFileName, append = TRUE)   
  return(result)
}
