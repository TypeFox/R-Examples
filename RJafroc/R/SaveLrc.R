#' @importFrom tools file_ext
SaveLrc <- function(dataset, fileName, dataDscrpt) {
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
  if (dataset$dataType != "ROC") {
    stop("Only ROC data file can be saved as *.lrc format.")
  }
  fileExt <- file_ext(fileName)
  if (fileExt != "lrc") {
    stop("The extension of LRC file name must be *.lrc")
  }
  write(dataDscrpt, fileName)
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
  for (j in 1:J) {
    write(readerID[j], fileName, append = TRUE)
    if (j == 1) {
      string <- ""
      for (i in 1:I) {
        string <- paste0(string, sprintf("   \"%s\"", modalityID[i]))
      }
      write(string, fileName, append = TRUE)
      
      string <- ""
      for (i in 1:I) {
        string <- paste0(string, "       L")
      }
      write(string, fileName, append = TRUE)
    }
    
    for (k in 1:K1) {
      string <- ""
      for (i in 1:I) {
        if (NL[i, j, k, 1] != UNINITIALIZED) {
          string <- paste0(string, sprintf("   %f", NL[i, j, k, 1]))
        } else {
          string <- paste0(string, sprintf("   %f", 0))
        }
      }
      string <- paste0(string, sprintf("   Non-diseased case %d", k))
      write(string, fileName, append = TRUE)
    }
    write("*", fileName, append = TRUE)
    
    for (k in 1:K2) {
      string <- ""
      for (i in 1:I) {
        if (LL[i, j, k, 1] != UNINITIALIZED) {
          string <- paste0(string, sprintf("   %f", LL[i, j, k, 1]))
        } else {
          string <- paste0(string, sprintf("   %f", 0))
        }
      }
      string <- paste0(string, sprintf("   Diseased Case %d", k))
      write(string, fileName, append = TRUE)
    }
    write("*", fileName, append = TRUE)
  }
  write("#", fileName, append = TRUE)
} 
