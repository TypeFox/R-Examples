mainNetFunction <-
function(counts, adjMat, nchips, plotPath = "", tfList = NULL){
  message("Evaluating Mutual Estimate Methods")
  miMLTime <- system.time(miML <- parMIEstimate(counts, method = "ML", unit = "nat", nchips = nchips, tfList = tfList))
  print(paste("ML Estimate Time:", round(miMLTime["elapsed"], 2)))
  miMMTime <- system.time(miMM <- parMIEstimate(counts, method = "MM", unit = "nat", nchips = nchips, tfList = tfList))
  print(paste("MM Estimate Time:", round(miMMTime["elapsed"], 2)))
  miBJTime <- system.time(miBJ <- parMIEstimate(counts, method = "Bayes", unit = "nat", nchips = nchips, priorHyperParam = "Jeffreys", tfList = tfList))
  print(paste("BJ Estimate Time:", round(miBJTime["elapsed"], 2)))
  miBBTime <- system.time(miBB <- parMIEstimate(counts, method = "Bayes", unit = "nat", nchips = nchips, priorHyperParam = "BLUnif", tfList = tfList))
  print(paste("BB Estimate Time:", round(miBBTime["elapsed"], 2)))
  miBPTime <- system.time(miBP <- parMIEstimate(counts, method = "Bayes", unit = "nat", nchips = nchips, priorHyperParam = "Perks", tfList = tfList))
  print(paste("BP Estimate Time:", round(miBPTime["elapsed"], 2)))
  miBMTime <- system.time(miBM <- parMIEstimate(counts, method = "Bayes", unit = "nat", nchips = nchips, priorHyperParam = "MiniMax", tfList = tfList))
  print(paste("BM Estimate Time:", round(miBMTime["elapsed"], 2)))
  miCSTime <- system.time(miCS <- parMIEstimate(counts, method = "CS", unit = "nat", nchips = nchips, tfList = tfList))
  print(paste("CS Estimate Time:", round(miCSTime["elapsed"], 2)))
  miSHTime <- system.time(miSH <- parMIEstimate(counts, method = "Shrink", unit = "nat", nchips = nchips, tfList = tfList))
  print(paste("SH Estimate Time:", round(miSHTime["elapsed"], 2)))
  miKDTime <- system.time(miKD <- parMIEstimate(counts, method = "KD", nchips = nchips, tfList = tfList))
  print(paste("KD Estimate Time:", round(miKDTime["elapsed"], 2)))
  miKNNTime <- system.time(miKNN <- parMIEstimate(counts, method = "KNN", unit = "nat", k = 3, nchips = nchips, tfList = tfList))
  print(paste("KNN Estimate Time:", round(miKNNTime["elapsed"], 2)))
  
  miEst <- list(miML = miML, miMM = miMM, miBJ = miBJ, miBB = miBB, miBP = miBP, miBM = miBM, miCS = miCS, miSH = miSH, miKD = miKD, miKNN = miKNN)
  
  message("Calculating Performance Indexes")
  valML <- performanceIndex(miML, adjMat)
  valMM <- performanceIndex(miMM, adjMat)
  valBJ <- performanceIndex(miBJ, adjMat)
  valBB <- performanceIndex(miBB, adjMat)
  valBP <- performanceIndex(miBP, adjMat)
  valBM <- performanceIndex(miBM, adjMat)
  valCS <- performanceIndex(miCS, adjMat)
  valSH <- performanceIndex(miSH, adjMat)
  valKD <- performanceIndex(miKD, adjMat)
  valKNN <- performanceIndex(miKNN, adjMat)
  
  valMet <- list(valML = valML, valMM = valMM, valBJ = valBJ, valBB = valBB, valBP = valBP, valBM = valBM, valCS = valCS, valSH = valSH, valKD = valKD, valKNN = valKNN)
  
  message("Generating ROC and PR Curves")
  png(paste(plotPath, "allROC.png", sep = ""), width = 855, height = 543)
  ccol <- c("red", "blue", "green", "black", "gray", "yellow", "aquamarine3", "darkmagenta", "burlywood4", "orange")
  plot(c(0, valML[, "FPR"], 1), c(0, valML[, "Recall"], 1), type = "l", col = ccol[1], xlab = "FP rate", ylab = "TP rate", main = "ROC Curve", xlim = 0:1, ylim = 0:1)
  lines(c(0, valMM[, "FPR"], 1), c(0, valMM[, "Recall"], 1), type = "l", col = ccol[2])
  lines(c(0, valBJ[, "FPR"], 1), c(0, valBJ[, "Recall"], 1), type = "l", col = ccol[3])
  lines(c(0, valBB[, "FPR"], 1), c(0, valBB[, "Recall"], 1), type = "l", col = ccol[4])
  lines(c(0, valBP[, "FPR"], 1), c(0, valBP[, "Recall"], 1), type = "l", col = ccol[5])
  lines(c(0, valBM[, "FPR"], 1), c(0, valBM[, "Recall"], 1), type = "l", col = ccol[6])
  lines(c(0, valCS[, "FPR"], 1), c(0, valCS[, "Recall"], 1), type = "l", col = ccol[7])
  lines(c(0, valSH[, "FPR"], 1), c(0, valSH[, "Recall"], 1), type = "l", col = ccol[8])
  lines(c(0, valKD[, "FPR"], 1), c(0, valKD[, "Recall"], 1), type = "l", col = ccol[9])
  lines(c(0, valKNN[, "FPR"], 1), c(0, valKNN[, "Recall"], 1), type = "l", col = ccol[10])
  lines(0:1, 0:1, col = "black")
  legend("bottomright", inset = .05, legend = c("ML", "MM", "BJ", "BB", "BP", "BM", "CS", "SH", "KD", "KNN"), title = "MI Estimators:", fill = ccol)
  dev.off()
  
  png(paste(plotPath, "allPR.png", sep = ""), width = 855, height = 543)
  ccol <- c("red", "blue", "green", "black", "gray", "yellow", "aquamarine3", "darkmagenta", "burlywood4", "orange")
  plot(c(0, valML[, "Recall"]), c(0, valML[, "Precision"]), type = "l", col = ccol[1], xlab = "recall", ylab = "precision", main = "PR Curve", xlim = 0:1, ylim = 0:1)
  lines(c(0, valMM[, "Recall"]), c(0, valMM[, "Precision"]), type = "l", col = ccol[2])
  lines(c(0, valBJ[, "Recall"]), c(0, valBJ[, "Precision"]), type = "l", col = ccol[3])
  lines(c(0, valBB[, "Recall"]), c(0, valBB[, "Precision"]), type = "l", col = ccol[4])
  lines(c(0, valBP[, "Recall"]), c(0, valBP[, "Precision"]), type = "l", col = ccol[5])
  lines(c(0, valBM[, "Recall"]), c(0, valBM[, "Precision"]), type = "l", col = ccol[6])
  lines(c(0, valCS[, "Recall"]), c(0, valCS[, "Precision"]), type = "l", col = ccol[7])
  lines(c(0, valSH[, "Recall"]), c(0, valSH[, "Precision"]), type = "l", col = ccol[8])
  lines(c(0, valKD[, "Recall"]), c(0, valKD[, "Precision"]), type = "l", col = ccol[9])
  lines(c(0, valKNN[, "Recall"]), c(0, valKNN[, "Precision"]), type = "l", col = ccol[10])
  legend("topright", inset = .05, legend = c("ML", "MM", "BJ", "BB", "BP", "BM", "CS", "SH", "KD", "KNN"), title = "MI Estimators:", fill = ccol)
  dev.off()
  
  message("Creating Result Table")
  resTable <- rbind(valML[which(valML$Fscore == max(valML$Fscore))[1], ],
                    valMM[which(valMM$Fscore == max(valMM$Fscore))[1], ],
                    valBJ[which(valBJ$Fscore == max(valBJ$Fscore))[1], ],
                    valBB[which(valBB$Fscore == max(valBB$Fscore))[1], ],
                    valBP[which(valBP$Fscore == max(valBP$Fscore))[1], ],
                    valBM[which(valBM$Fscore == max(valBM$Fscore))[1], ],
                    valCS[which(valCS$Fscore == max(valCS$Fscore))[1], ],
                    valSH[which(valSH$Fscore == max(valSH$Fscore))[1], ],
                    valKNN[which(valKD$Fscore == max(valKD$Fscore))[1], ],
                    valKNN[which(valKNN$Fscore == max(valKNN$Fscore))[1], ]) #check when more than one max Fscore!!!
  rownames(resTable) <- c("ML", "MM", "BJ", "BB", "BP", "BM", "CS", "SH", "KD", "KNN")
  
  AUROC <- c(aucDisc(valML[, "FPR"], valML[, "Recall"]),
             aucDisc(valMM[, "FPR"], valMM[, "Recall"]),
             aucDisc(valBJ[, "FPR"], valBJ[, "Recall"]),
             aucDisc(valBB[, "FPR"], valBB[, "Recall"]),
             aucDisc(valBP[, "FPR"], valBP[, "Recall"]),
             aucDisc(valBM[, "FPR"], valBM[, "Recall"]),
             aucDisc(valCS[, "FPR"], valCS[, "Recall"]),
             aucDisc(valSH[, "FPR"], valSH[, "Recall"]),
             aucDisc(valKD[, "FPR"], valKD[, "Recall"]),
             aucDisc(valKNN[, "FPR"], valKNN[, "Recall"]))
  AUPR <- c(aucDisc(valML[, "Recall"], valML[, "Precision"]),
            aucDisc(valMM[, "Recall"], valMM[, "Precision"]),
            aucDisc(valBJ[, "Recall"], valBJ[, "Precision"]),
            aucDisc(valBB[, "Recall"], valBB[, "Precision"]),
            aucDisc(valBP[, "Recall"], valBP[, "Precision"]),
            aucDisc(valBM[, "Recall"], valBM[, "Precision"]),
            aucDisc(valCS[, "Recall"], valCS[, "Precision"]),
            aucDisc(valSH[, "Recall"], valSH[, "Precision"]),
            aucDisc(valKD[, "Recall"], valKD[, "Precision"]),
            aucDisc(valKNN[, "Recall"], valKNN[, "Precision"]))
  
  resTable <- cbind(resTable, AUROC, AUPR)
  
  ans <- list(miEst = miEst, valMet = valMet, resTable = resTable)
  return(ans)
}
