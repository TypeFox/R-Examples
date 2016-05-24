#' @importFrom tools file_ext
#' @importFrom stringr str_trim
ReadImrmc <- function(fileName) {
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
  fileExt <- file_ext(fileName)
  if (fileExt != "imrmc") {
    stop("The extension of iMRMC file name must be *.imrmc")
  }
  fileLines <- str_trim(readLines(fileName))
  lineNum <- 1
  
  for (n in 1:4) {
    lineTxt <- unlist(strsplit(fileLines[lineNum], split = ":", fixed = TRUE))
    while (!(lineTxt[1] %in% c("N0", "N1", "NR", "NM"))) {
      lineNum <- lineNum + 1
      lineTxt <- unlist(strsplit(fileLines[lineNum], split = ":", fixed = TRUE))
    }
    if (lineTxt[1] == "N0") {
      K1 <- as.numeric(lineTxt[2])
    } else if (lineTxt[1] == "N1") {
      K2 <- as.numeric(lineTxt[2])
    } else if (lineTxt[1] == "NM") {
      I <- as.numeric(lineTxt[2])
    } else if (lineTxt[1] == "NR") {
      J <- as.numeric(lineTxt[2])
    }
    lineNum <- lineNum + 1
  }
  K <- K1 + K2
  
  while (fileLines[lineNum] != "BEGIN DATA:") lineNum <- lineNum + 1
  lineNum <- lineNum + 1
  while (fileLines[lineNum] == "") lineNum <- lineNum + 1
  dataTable <- NULL
  while (length(lineTxt <- unlist(strsplit(fileLines[lineNum], split = ",", fixed = TRUE))) == 4) {
    dataTable <- rbind(dataTable, lineTxt)
    lineNum <- lineNum + 1
  }
  
  caseId <- unique(as.numeric(dataTable[, 2]))
  if (any(is.na(caseId))) {
    errMsg <- "Case IDs must be integers."
    stop(errMsg)
  }
  if (any(is.na(as.numeric(dataTable[, 4])))) {
    errMsg <- "Ratings must be numeric."
    stop(errMsg)
  }
  if (length(caseId) != K) {
    errMsg <- "The number of cases in the dataset is different from the number in the study description."
    stop(errMsg)
  }
  
  truthTable <- dataTable[as.numeric(dataTable[, 1]) == -1, ]
  if (any(!(caseId %in% as.numeric(truthTable[, 2])))) {
    errMsg <- sprintf("Case %d cannot be found in rows recording truth states.", caseId[!(caseId %in% as.numeric(truthTable[, 2]))])
    stop(errMsg)
  }
  
  if (any(!(as.numeric(truthTable[, 4]) %in% c(0, 1)))) {
    errMsg <- "Cases' truth states must be 0 or 1."
    stop(errMsg)
  }
  normalCases <- unique(as.numeric(truthTable[as.numeric(truthTable[, 4]) == 0, 2]))
  abnormalCases <- unique(as.numeric(truthTable[as.numeric(truthTable[, 4]) == 1, 2]))
  
  fpTable <- dataTable[(as.numeric(dataTable[, 2]) %in% normalCases) & as.numeric(dataTable[, 1]) != -1, ]
  tpTable <- dataTable[(as.numeric(dataTable[, 2]) %in% abnormalCases) & as.numeric(dataTable[, 1]) != -1, ]
  
  readerID <- as.character(sort(unique(c(fpTable[, 1], tpTable[, 1]))))
  if (J != length(readerID)) {
    errMsg <- "The number of readers in the dataset is different from the number in the study description."
    stop(errMsg)
  }
  
  modalityID <- as.character(sort(unique(c(fpTable[, 3], tpTable[, 3]))))
  if (I != length(modalityID)) {
    errMsg <- "The number of modalities in the dataset is different from the number in the study description."
    stop(errMsg)
  }
  
  NL <- array(UNINITIALIZED, dim = c(I, J, K, 1))
  LL <- array(UNINITIALIZED, dim = c(I, J, K2, 1))
  for (i in 1:I) {
    for (j in 1:J) {
      for (k1 in 1:K1) {
        caseIndx <- which((fpTable[, 1] == readerID[j]) & (fpTable[, 3] == modalityID[i]) & (as.numeric(fpTable[, 2]) == normalCases[k1]))
        NL[i, j, k1, 1] <- as.numeric(fpTable[caseIndx, 4])
      }
      for (k2 in 1:K2) {
        caseIndx <- which((tpTable[, 1] == readerID[j]) & (tpTable[, 3] == modalityID[i]) & (as.numeric(tpTable[, 2]) == abnormalCases[k2]))
        LL[i, j, k2, 1] <- as.numeric(tpTable[caseIndx, 4])
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
