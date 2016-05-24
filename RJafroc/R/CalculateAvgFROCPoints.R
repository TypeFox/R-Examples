CalculateAvgFROCPoints <- function(NL, LL, modalityID, readerID, maxNL, lesionNum, plottingModalities, plottingReaders) {
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
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
  
  sumLL <- sum(lesionNum)
  
  NLF <- list(NULL)
  LLF <- list(NULL)
  avgIndex <- 1
  for (i in 1:I) {
    for (j in 1:J) {
      nl <- NL[i, j, , ][NL[i, j, , ] != UNINITIALIZED]
      ll <- LL[i, j, , ][LL[i, j, , ] != UNINITIALIZED]
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
      
      nlf <- cumsum(as.vector(nlTable))/K
      llf <- cumsum(as.vector(llTable))/sumLL
      
      NLF[[avgIndex]] <- 0:length(scores)
      LLF[[avgIndex]] <- NLF[[avgIndex]]
      
      numNL <- 1
      numLL <- 1
      for (k in 1:length(scores)) {
        if (!is.na(nlScores[numNL]) && nlScores[numNL] >= scores[k]) {
          NLF[[avgIndex]][k + 1] <- nlf[numNL]
          numNL <- numNL + 1
        } else {
          NLF[[avgIndex]][k + 1] <- NLF[[avgIndex]][k]
        }
        
        if (!is.na(llScores[numLL]) && llScores[numLL] >= scores[k]) {
          LLF[[avgIndex]][k + 1] <- llf[numLL]
          numLL <- numLL + 1
        } else {
          LLF[[avgIndex]][k + 1] <- LLF[[avgIndex]][k]
        }
      }
      avgIndex <- avgIndex + 1
    }
  }
  
  maxNLF <- max(sapply(NLF, max))
  
  plotAvgStep <- 0.005
  avgNLF <- seq(0, (1 + plotAvgStep), plotAvgStep) * maxNLF
  avgLLFArray <- array(dim = c(I * J, length(avgNLF)))
  avgIndex <- 1
  for (i in 1:I) {
    for (j in 1:J) {
      tempLLFAvg <- avgNLF
      oldLLFIndx <- 1
      NLF[[avgIndex]] <- c(NLF[[avgIndex]], maxNLF)
      LLF[[avgIndex]] <- c(LLF[[avgIndex]], LLF[[avgIndex]][length(LLF[[avgIndex]])])
      for (k in 2:length(LLF[[avgIndex]])) {
        LLFIndx <- rev(which(NLF[[avgIndex]][k] >= avgNLF))[1] + 1
        if (!is.na(tempLLFAvg[LLFIndx])) {
          tempLLFAvg[LLFIndx] <- LLF[[avgIndex]][k]
          if ((LLFIndx - 1) >= (oldLLFIndx + 1)) {
            tempLLFAvg[(oldLLFIndx + 1):(LLFIndx - 1)] <- tempLLFAvg[oldLLFIndx] + seq(from = plotAvgStep, by = plotAvgStep, length.out = (LLFIndx - oldLLFIndx - 1)) * maxNLF * ((tempLLFAvg[(LLFIndx)] - 
                                                                                                                                                                                     tempLLFAvg[(oldLLFIndx)])/((LLFIndx - oldLLFIndx - 1) * plotAvgStep * maxNLF))
          }
          oldLLFIndx <- LLFIndx
        } else {
          tempLLFAvg[LLFIndx - 2] <- LLF[[avgIndex]][k]
        }
      }
      
      if (LLFIndx < 1/plotAvgStep) {
        tempLLFAvg[(LLFIndx + 1):((1/plotAvgStep) + 1)] <- NA
      }
      
      avgLLFArray[avgIndex, ] <- tempLLFAvg
      avgIndex <- avgIndex + 1
    }
  }
  
  LLFAvg <- rep(NA, length(avgNLF))
  for (p in 1:length(avgNLF)) LLFAvg[p] <- mean(avgLLFArray[!is.na(avgLLFArray[, p]), p])
  LLFAvg <- LLFAvg[1:sum(!is.na(LLFAvg))]
  avgNLF <- avgNLF[1:length(LLFAvg)]
  
  class <- paste(paste("M-"), paste(modalityID, collapse = " "), "\n", paste("R-"), paste(readerID, collapse = " "), sep = "")
  FROCPoints <- data.frame(NLF = avgNLF, LLF = LLFAvg, class = class, type = "averaged")
  return(FROCPoints)
} 
