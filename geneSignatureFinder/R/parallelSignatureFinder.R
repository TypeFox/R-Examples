parallelSignatureFinder <-
function(cpuCluster, startingGene, 
	  logFileName = "",
	  coeffMissingAllowed = 0.75,
    subsetToUse = 1:ncol(geData), zero = NULL) {

  mmessage = FALSE
  n <- nrow(geData)
  m <- ncol(geData)

  toExplore <- rep(FALSE, m)
  toExplore[subsetToUse]  <- TRUE
  toExplore[startingGene]  <- FALSE
  
#  logFileName <- paste(logFileName, "SignatureFinderLog.txt", sep = "")
#  cat(paste("Working on ", sum(toExplore), " genes observed on ", n, " samples.\n", sep = ""), file = logFileName, append = FALSE)
#  cat(paste("Starting on ", ttime <- Sys.time(), "\n", sep = ""), file = logFileName, append = TRUE)
#  cat(paste("Starting signature: ", paste(colnames(geData)[startingGene], collapse = ", "), ";\n", sep = ""), file = logFileName, append = TRUE)
  
  
  message(paste("\n\nWorking on ", sum(toExplore), " genes observed on ", n, " samples.", sep = ""))
  message(paste("Starting on ", ttime <- Sys.time(), sep = ""))
  message(paste("Starting signature: ", paste(colnames(geData)[startingGene], collapse = ", " ), ";", sep = ""))
  
  aClassify <- rep(NA, n)
  ssignature <- startingGene
  
  runs <- 0 
  result <- list()  
  result$signatureName <-  (colnames(geData)[startingGene])[1]
  result$startingSignature <-  colnames(geData)[startingGene]
  result$coeffMissingAllowed <- coeffMissingAllowed
  ####### aggiunto il 14/10/2013
  result$zero <- zero

######## da rivedere; si attiva nel caso lo startingGene sia una sequenza
  if(length(ssignature) > 1) 
    notMissing <- apply(!is.na(geData[, ssignature]), 1, sum) else
  notMissing <- !is.na(geData[, ssignature])
  notMissing <- notMissing > 0


  #clusters <- classify(geData[notMissing, ssignature])$clusters  

  clusters <- classify(geData[, ssignature])$clusters  
  ####### modificato il 14/10/2013; qua si verifichera' un errore se ci sono missing
  ####### perche' stData non risultera' allineata con i missing
  clusters <- goodAndPoorClassification(clusters)
  
  runningDistance <- survdiff(stData[notMissing] ~ clusters[notMissing])$chisq

  result$startingTValue <- runningDistance
  result$startingPValue <- 1 - pchisq(runningDistance, df = 1)  
  ####### modificato il 14/10/2013;
  #tmpClassification <- rep(NA, n)
  #tmpClassification[notMissing] <- clusters 
  #result$startingClassification <- as.factor(tmpClassification)
  #levels(result$startingClassification) <- c("good", "poor")
  result$startingClassification <- clusters

  #cat(paste("tValue = ", runningDistance, " (", result$startingPValue, ")\n", sep = ""), file = logFileName, append = TRUE)
  message(paste("tValue = ", runningDistance, " (", result$startingPValue, ")\n", sep = ""))
####################### end: da rivedere

  
  tmpTime <- Sys.time()
  exitFromMain <- FALSE
  repeat {#### MAIN LOOP #######################
    if(mmessage) message(runs+1, " MAIN LOOP: ", exitFromMain)
    if(exitFromMain) break
    runs <- runs + 1
    distances <- rep(0, m) #variato
    nnrow = sum(toExplore)###<-
    signatureToExplore <- matrix(ssignature,
                                 ncol = length(ssignature),
                                 nrow = nnrow,
                                 byrow = TRUE)    
    signatureToExplore <- cbind(which(toExplore), signatureToExplore)
    
    if(nnrow == 1) {###<-
      distances[toExplore] <- tValueFun(signatureToExplore) ###<-
      exitFromMain <- TRUE###<-
    } else###<-
      distances[toExplore] <- parApply(cpuCluster, signatureToExplore, 1, tValueFun, ###<-
                                     coeffMissingAllowed = coeffMissingAllowed)###<-
      
    
    #howManyCandidates <- sum(distances > runningDistance)
    #message("sum(distances > runningDistance) = ", howManyCandidates)
    #foundNewGene  <- FALSE
    exitFromInner <-  FALSE
    repeat {### INNER LOOP ##########################
      if(mmessage) message("MAIN = ", runs, ", INNER LOOP: ", exitFromInner, " time: ", round(Sys.time() - tmpTime, 2))
      if(exitFromInner) break
      maxDistance <- max(distances) #variato
      if(maxDistance >= runningDistance) { #variato
        if(mmessage) message("if(maxDistance >= runningDistance), time: ", round(Sys.time() - tmpTime, 2))
        candidate <- which(distances == maxDistance)

        if(length(candidate) == 1) {
          if(mmessage) message(" if(length(candidate) == 1), time: ", round(Sys.time() - tmpTime, 2), " ", candidate)
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
          #############################
          #################################################################################
          if(is.null(zero)) checks <- checkOne & checkTwo else {
            mean1 <- mean(geData[clusters == 1, candidate])
            mean2 <- mean(geData[clusters == 2, candidate])
            checkThree  <- (mean1- zero) * (mean2- zero) < 0
            checks <- checkOne & checkTwo & checkThree                
          }
          if(checks) {
            if(mmessage) message("if(checks), time: ", round(Sys.time() - tmpTime, 2))
            ssignature <- c(ssignature, candidate)
            toExplore[ssignature] <- FALSE
            runningDistance <- maxDistance
            tmp <- paste("... improved signature: ",
                         paste(colnames(geData)[ssignature], collapse = ", "),
                         ";\ntValue = ", runningDistance,  " (",
                         1 - pchisq(runningDistance, df = 1), ")\n", sep = "")
            message(tmp)
            
            exitFromInner <- TRUE
          } # end of if(checkOne & checkTwo & checkThree) {
        } else { #   if(length(candidate) == 1)
          if(mmessage) message(" if(length(candidate) == 1) else, time: ", round(Sys.time() - tmpTime, 2))
          if(length(candidate) > 0.01*m) {
            if(mmessage) message(" if(length(candidate) > 0.01*m)")
            exitFromInner <- TRUE
            exitFromMain <- TRUE
            break
          }
          tmpCandidate <- parallelFindBest(cpuCluster, runningDistance, ssignature, candidate)
          if(length(tmpCandidate) > 0) { 
            if(mmessage) message("if(length(tmpCandidate) > 0), ", length(tmpCandidate), ",  time: ", round(Sys.time() - tmpTime, 2))
            tmp <- survdiff(stData ~ classify(geData[, c(ssignature, tmpCandidate)])$clusters)$chisq
            if(tmp > runningDistance) {
              if(mmessage) message(" if(tmp > runningDistance), time: ", round(Sys.time() - tmpTime, 2))
              # controlli inseriti il 14/10/2013
              ################################################################################
              clusters <- classify(geData[notMissing, c(ssignature, tmpCandidate)])$clusters
              checkOne <- (min(table(clusters)) > floor(0.1 * n))
              # qua si controlla che le curve di sopravvivenza non si incrocino
              sf0 <- survfit(stData[clusters == 1] ~ 1)$surv
              sf1 <- survfit(stData[clusters == 2] ~ 1)$surv
              checkTwo <- sum(fivenum(sf0) > fivenum(sf1))
              checkTwo <- ((checkTwo == 0) | (checkTwo == 5))
              #############################
              #################################################################################
              if(is.null(zero)) checks <- checkOne & checkTwo else {
                mean1 <- mean(geData[clusters == 1, tmpCandidate])
                mean2 <- mean(geData[clusters == 2, tmpCandidate])
                checkThree  <- (mean1- zero) * (mean2- zero) < 0
                checks <- checkOne & checkTwo & checkThree                
              }
              if(checks) {
                if(mmessage) message("if(checks), time: ", round(Sys.time() - tmpTime, 2))
              #################################
                candidate <- tmpCandidate
                ssignature <- c(ssignature, candidate)
                runningDistance <- tmp
                toExplore[ssignature] <- FALSE
                tmp <- paste("... improved signature: ",
                             paste(colnames(geData)[ssignature], collapse = ", "),
                             ";\ntValue = ", runningDistance,  " (",
                             1 - pchisq(runningDistance, df = 1), ")\n", sep = "")
                message(tmp)
                exitFromInner <- TRUE
              } else {
                if(mmessage) message("if(checks) else, time: ", round(Sys.time() - tmpTime, 2))
                #toExplore[tmpCandidate] <- FALSE   
                distances[tmpCandidate] <- 0 #variato                     
              }
              #######################
            }  # end of  if(tmp > runningDistance)
          } # end of if(length(tmpCandidate) > 0)      
        } # end of if(length(candidate) == 1) else ...
        distances[candidate] <- 0 #variato        
      } else { #  if(maxDistance >= runningDistance)
        exitFromInner <- TRUE
        exitFromMain <- TRUE
      } # end of  if(maxDistance >= runningDistance) else
    } ### END of INNER LOOP ###########
          if(mmessage) message("END of INNER LOOP, time: ", round(Sys.time() - tmpTime, 2))
  } ### END of MAIN LOOP ###########
  if(mmessage) message("END of MAIN LOOP, time: ", round(Sys.time() - tmpTime, 2))
  
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
  #cat(paste("\n\nfinal signature: ", paste(colnames(geData)[ssignature], collapse = " "), sep = ""), file = logFileName, append = TRUE)
  #cat(paste("\ntValue = ", runningDistance, " (", 1 - pchisq(runningDistance, df = 1), ")\n",  sep = ""), file = logFileName, append = TRUE)
  #cat(paste("\nlength of the signature = ", length(ssignature), sep = ""), file = logFileName, append = TRUE)
  #cat(paste("\nnumber of joint missing values = ", sum(!notMissing), " (", 100*round(sum(!notMissing)/n,2), "%)", sep = ""), file = logFileName, append = TRUE)
  #cat(paste("\n\nEnd of computation at ", t2 <- Sys.time(), "; elapsed time: ", t2 - ttime, ".", sep = ""), file = logFileName, append = TRUE) 
  
  
  message(paste("final signature: ", paste(colnames(geData)[ssignature], collapse = " "), sep = ""))
  message(paste("tValue = ", runningDistance, " (", 1 - pchisq(runningDistance, df = 1), ")\n",  sep = ""))
  message(paste("length of the signature = ", length(ssignature), sep = ""))
  message(paste("number of joint missing values = ", sum(!notMissing), " (", 100*round(sum(!notMissing)/n,2), "%)", sep = ""))
  message(paste("End of computation at ", t2 <- Sys.time(), "; elapsed time: ", t2 - ttime, ".\n\n", sep = ""))   
  return(result)
}
