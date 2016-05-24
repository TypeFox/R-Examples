#' @importFrom tools file_ext
SaveOrDbmMrmc <- function(dataset, fileName) {
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
  if (dataset$dataType != "ROC") {
    stop("Only ROC data file can be saved in MRMC format.")
  }
  fileExt <- file_ext(fileName)
  if (!fileExt %in% c("csv", "lrc", "txt")) {
    stop("The extension of MRMC file name must be *.txt or *.csv or *.lrc.")
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
  
  write("reader,treatment,case,truth,rating", fileName)
  for (i in 1:I) {
    for (j in 1:J) {
      for (k in 1:K1) {
        if (NL[i, j, k, 1] != UNINITIALIZED) {
          write(sprintf("%s,%s,%d,%d,%f", dataset$readerID[j], dataset$modalityID[i], k, 0, NL[i, j, k, 1]), fileName, append = TRUE)
        } else {
          write(sprintf("%s,%s,%d,%d,%f", dataset$readerID[j], dataset$modalityID[i], k, 0, -2000), fileName, append = TRUE)
        }
      }
      
      for (k in 1:K2) {
        if (LL[i, j, k, 1] != UNINITIALIZED) {
          write(sprintf("%s,%s,%d,%d,%f", dataset$readerID[j], dataset$modalityID[i], k + K1, 1, LL[i, j, k, 1]), fileName, append = TRUE)
        } else {
          write(sprintf("%s,%s,%d,%d,%f", dataset$readerID[j], dataset$modalityID[i], k + K1, 1, -2000), fileName, append = TRUE)
        }
      }
    }
  }
} 
