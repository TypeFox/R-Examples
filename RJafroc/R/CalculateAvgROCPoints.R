CalculateAvgROCPoints <- function(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities, plottingReaders) {
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
  
  plotAvgStep <- 0.005
  FPFAvg <- seq(0, 1, plotAvgStep)
  avgTPF <- rep(0, length(FPFAvg))
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
      
      
      tempAvgTPF <- FPFAvg
      oldTPFindx <- 1
      for (k in 2:length(scores)) {
        TPFindx <- rev(which(FPF[k] >= FPFAvg))[1] + 1
        if (!is.na(tempAvgTPF[TPFindx])) {
          tempAvgTPF[TPFindx] <- TPF[k]
          if ((TPFindx - 1) >= (oldTPFindx + 1)) {
            tempAvgTPF[(oldTPFindx + 1):(TPFindx - 1)] <- tempAvgTPF[oldTPFindx] + seq(from = plotAvgStep, by = plotAvgStep, length.out = (TPFindx - oldTPFindx - 1)) * ((tempAvgTPF[(TPFindx)] - tempAvgTPF[(oldTPFindx)])/((TPFindx - 
                                                                                                                                                                                                                                oldTPFindx - 1) * plotAvgStep))
          }
          oldTPFindx <- TPFindx
        } else {
          tempAvgTPF[TPFindx - 2] <- TPF[k]
        }
      }
      if (TPFindx < 1/plotAvgStep) {
        oldTPFindx <- TPFindx
        TPFindx <- (1/plotAvgStep) + 1
        tempAvgTPF[(oldTPFindx + 1):(TPFindx - 1)] <- tempAvgTPF[oldTPFindx] + seq(from = plotAvgStep, by = plotAvgStep, length.out = (TPFindx - oldTPFindx - 1)) * ((tempAvgTPF[(TPFindx)] - tempAvgTPF[(oldTPFindx)])/((TPFindx - 
                                                                                                                                                                                                                            oldTPFindx - 1) * plotAvgStep))
      }
      
      avgTPF <- avgTPF + tempAvgTPF
    }
  }
  avgTPF <- avgTPF/(I * J)
  
  class <- paste(paste("M-"), paste(modalityID, collapse = " "), "\n", paste("R-"), paste(readerID, collapse = " "), sep = "")
  ROCPoints <- data.frame(FPF = FPFAvg, TPF = avgTPF, class = class, type = "averaged")
  return(ROCPoints)
} 
