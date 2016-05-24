#' @importFrom tools file_ext
SaveImrmc <- function(dataset, fileName, dataDscrpt) {
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
  if (dataset$dataType != "ROC") 
    stop("Only ROC data file can be saved as iMRMC format.")
  fileExt <- file_ext(fileName)
  if (fileExt != "imrmc") {
    stop("The extension of iMRMC file name must be *.imrmc.")
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
  write(dataDscrpt, fileName)
  write(sprintf("N0:%d", K1), fileName, append = TRUE)
  write(sprintf("N1:%d", K2), fileName, append = TRUE)
  write(sprintf("NR:%d", J), fileName, append = TRUE)
  write(sprintf("NM:%d", I), fileName, append = TRUE)
  write(sprintf("BEGIN DATA:"), fileName, append = TRUE)
  for (k in 1:K1) {
    write(sprintf("-1,%d,truth,0", k), fileName, append = TRUE)
  }
  for (k in 1:K2) {
    write(sprintf("-1,%d,truth,1", k + K1), fileName, append = TRUE)
  }
  for (i in 1:I) {
    for (j in 1:J) {
      for (k in 1:K1) {
        if (NL[i, j, k, 1] != UNINITIALIZED) {
          write(sprintf("%s,%d,%s,%f", dataset$readerID[j], k, dataset$modalityID[i], NL[i, j, k, 1]), fileName, append = TRUE)
        } else {
          write(sprintf("%s,%d,%s,%f", dataset$readerID[j], k, dataset$modalityID[i], -2000), fileName, append = TRUE)
        }
      }
      
      for (k in 1:K2) {
        if (LL[i, j, k, 1] != UNINITIALIZED) {
          write(sprintf("%s,%d,%s,%f", dataset$readerID[j], k + K1, dataset$modalityID[i], LL[i, j, k, 1]), fileName, append = TRUE)
        } else {
          write(sprintf("%s,%d,%s,%f", dataset$readerID[j], k + K1, dataset$modalityID[i], -2000), fileName, append = TRUE)
        }
      }
    }
  }
} 
