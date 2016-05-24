CalculateROCPoints <- function(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities, plottingReaders) {
  I <- length(plottingModalities)
  J <- dim(NL)[2]
  K <- dim(NL)[3]
  K2 <- dim(LL)[3]
  K1 <- K - K2
  
  if (all(is.character(plottingModalities))) 
    plottingModalities <- which(modalityID %in% plottingModalities)
  NL <- NL[plottingModalities, , , ]
  LL <- LL[plottingModalities, , , ]
  dim(NL) <- c(I, J, K, maxNL)
  dim(LL) <- c(I, J, K2, max(lesionNum))
  modalityID <- modalityID[plottingModalities]
  
  
  J <- length(plottingReaders)
  
  if (all(is.character(plottingReaders))) 
    plottingReaders <- which(readerID %in% plottingReaders)
  NL <- NL[, plottingReaders, , ]
  LL <- LL[, plottingReaders, , ]
  dim(NL) <- c(I, J, K, maxNL)
  dim(LL) <- c(I, J, K2, max(lesionNum))
  readerID <- readerID[plottingReaders]
  
  
  NL <- apply(NL, c(1, 2, 3), max)
  LL <- apply(LL, c(1, 2, 3), max)
  
  ROCPoints <- data.frame(FPF = NULL, TPF = NULL)
  for (i in 1:I) {
    for (j in 1:J) {
      nl <- NL[i, j, 1:K1]
      ll <- apply(cbind(NL[i, j, (K1 + 1):K], LL[i, j, ]), 1, max)
      
      nlTable <- rev(table(nl))
      llTable <- rev(table(ll))
      
      if (length(nlTable) == 1) {
        nlScores <- as.numeric(attr(nlTable, "names"))
      } else {
        nlScores <- as.numeric(unlist(attr(nlTable, "dimnames")))
      }
      
      if (length(llTable) == 1) {
        llScores <- as.numeric(attr(llTable, "names"))
      } else {
        llScores <- as.numeric(unlist(attr(llTable, "dimnames")))
      }
      scores <- sort(unique(c(nlScores, llScores)), decreasing = TRUE)
      
      fpf <- cumsum(as.vector(nlTable))/K1
      tpf <- cumsum(as.vector(llTable))/K2
      
      FPF <- 0:length(scores)
      TPF <- FPF
      
      numNL <- 1
      numLL <- 1
      for (k in 1:length(scores)) {
        if (!is.na(nlScores[numNL]) && nlScores[numNL] >= scores[k]) {
          FPF[k + 1] <- fpf[numNL]
          numNL <- numNL + 1
        } else {
          FPF[k + 1] <- FPF[k]
        }
        
        if (!is.na(llScores[numLL]) && llScores[numLL] >= scores[k]) {
          TPF[k + 1] <- tpf[numLL]
          numLL <- numLL + 1
        } else {
          TPF[k + 1] <- TPF[k]
        }
      }
      FPF[k + 1] <- 1
      TPF[k + 1] <- 1
      ROCPoints <- rbind(ROCPoints, data.frame(FPF = FPF, TPF = TPF, Modality = i, Reader = j))
    }
  }
  class <- paste("M-", modalityID[ROCPoints$Modality], "\n", "R-", readerID[ROCPoints$Reader], sep = "")
  ROCPoints <- data.frame(FPF = ROCPoints$FPF, TPF = ROCPoints$TPF, class = class, type = "individual")
  return(ROCPoints)
} 
