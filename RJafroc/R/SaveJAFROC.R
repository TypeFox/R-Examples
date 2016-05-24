#' @importFrom tools file_ext
SaveJAFROC <- function(dataset, fileName) {
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
  fileExt <- file_ext(fileName)
  if (!fileExt %in% c("xls", "xlsx")) {
    stop("The extension of JAFROC file name must be *.xls or *.xlsx.")
  }
  NL <- dataset$NL
  LL <- dataset$LL
  lesionNum <- dataset$lesionNum
  lesionID <- dataset$lesionID
  lesionWeight <- dataset$lesionWeight
  maxNL <- dim(NL)[4]
  dataType <- dataset$dataType
  modalityID <- dataset$modalityID
  readerID <- dataset$readerID
  I <- length(modalityID)
  J <- length(readerID)
  K <- dim(NL)[3]
  K2 <- dim(LL)[3]
  K1 <- K - K2
  
  caseIDs <- c(1:K1, rep(K1 + 1:K2, lesionNum))
  lesionIDs <- as.vector(t(lesionID))
  lesionIDs <- lesionIDs[lesionIDs != UNINITIALIZED]
  lesionIDs <- c(rep(0, K1), lesionIDs)
  lesionWeights <- as.vector(t(lesionWeight))
  lesionWeights <- lesionWeights[lesionWeights != UNINITIALIZED]
  lesionWeights <- c(rep(0, K1), lesionWeights)
  dataSheet <- data.frame(CaseID = as.integer(caseIDs), LesionID = as.integer(lesionIDs), Weight = lesionWeights)
  write.xlsx2(x = dataSheet, file = fileName, sheetName = "TRUTH", row.names = FALSE)
  
  dataSheet <- NULL
  for (i in 1:I) {
    for (j in 1:J) {
      for (k in 1:K) {
        for (l in 1:maxNL) {
          if (NL[i, j, k, l] != UNINITIALIZED) {
            dataSheet <- rbind(dataSheet, c(j, i, k, NL[i, j, k, l]))
          }
        }
      }
    }
  }
  dataSheet <- data.frame(ReaderID = readerID[dataSheet[, 1]], ModalityID = modalityID[dataSheet[, 2]], CaseID = as.integer(dataSheet[, 3]), NL_Rating = signif(dataSheet[, 4], 6))
  write.xlsx2(x = dataSheet, file = fileName, sheetName = "FP", row.names = FALSE, append = TRUE)
  
  dataSheet <- NULL
  for (i in 1:I) {
    for (j in 1:J) {
      for (k in 1:K2) {
        for (l in 1:lesionNum[k]) {
          if (LL[i, j, k, l] != UNINITIALIZED) {
            dataSheet <- rbind(dataSheet, c(j, i, k + K1, lesionID[k, l], LL[i, j, k, l]))
          }
        }
      }
    }
  }
  dataSheet <- data.frame(ReaderID = readerID[dataSheet[, 1]], ModalityID = modalityID[dataSheet[, 2]], CaseID = as.integer(dataSheet[, 3]), LesionID = as.integer(dataSheet[, 4]), LL_Rating = signif(dataSheet[, 5], 6))
  write.xlsx2(x = dataSheet, file = fileName, sheetName = "TP", row.names = FALSE, append = TRUE)
} 
