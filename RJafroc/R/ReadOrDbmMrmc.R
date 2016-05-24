#' @importFrom utils read.delim
ReadOrDbmMrmc <- function(fileName, delimiter) {
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
  dataTableFrame <- read.delim(fileName, sep = delimiter)
  dataTable <- NULL
  for (n in 1:5) {
    dataTable <- cbind(dataTable, as.character(dataTableFrame[, n]))
  }
  
  caseId <- unique(as.numeric(dataTable[, 3]))
  if (any(is.na(caseId))) {
    errMsg <- "Case IDs must be integers."
    stop(errMsg)
  }
  if (any(is.na(as.numeric(dataTable[, 5])))) {
    errMsg <- "Ratings must be integers."
    stop(errMsg)
  }
  if (any(!(as.numeric(dataTable[, 4]) %in% c(0, 1)))) {
    errMsg <- "Cases' truth states must be 0 or 1."
    stop(errMsg)
  }
  
  normalCases <- unique(as.numeric(dataTable[as.numeric(dataTable[, 4]) == 0, 3]))
  abnormalCases <- unique(as.numeric(dataTable[as.numeric(dataTable[, 4]) == 1, 3]))
  K1 <- length(normalCases)
  K2 <- length(abnormalCases)
  K <- K1 + K2
  
  fpTable <- dataTable[(as.numeric(dataTable[, 3]) %in% normalCases), ]
  tpTable <- dataTable[(as.numeric(dataTable[, 3]) %in% abnormalCases), ]
  
  readerID <- as.character(sort(unique(c(fpTable[, 1], tpTable[, 1]))))
  J <- length(readerID)
  
  modalityID <- as.character(sort(unique(c(fpTable[, 2], tpTable[, 2]))))
  I <- length(modalityID)
  
  NL <- array(UNINITIALIZED, dim = c(I, J, K, 1))
  LL <- array(UNINITIALIZED, dim = c(I, J, K2, 1))
  for (i in 1:I) {
    for (j in 1:J) {
      for (k1 in 1:K1) {
        caseIndx <- which((fpTable[, 1] == readerID[j]) & (fpTable[, 2] == modalityID[i]) & (as.numeric(fpTable[, 3]) == normalCases[k1]))
        NL[i, j, k1, 1] <- as.numeric(fpTable[caseIndx, 5])
      }
      for (k2 in 1:K2) {
        caseIndx <- which((tpTable[, 1] == readerID[j]) & (tpTable[, 2] == modalityID[i]) & (as.numeric(tpTable[, 3]) == abnormalCases[k2]))
        LL[i, j, k2, 1] <- as.numeric(tpTable[caseIndx, 5])
      }
    }
  }
  
  lesionNum <- rep(1, K2)
  lesionID <- array(1, dim = c(K2, 1))
  lesionWeight <- lesionID
  maxNL <- 1
  dataType <- "ROC"
  
  return(list(NL = NL, LL = LL, lesionNum = lesionNum, lesionID = lesionID, lesionWeight = lesionWeight, dataType = dataType, modalityID = modalityID, readerID = readerID))
} 
