#' @import xlsx
ReadJAFROC <- function(fileName) {
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
  wb <- loadWorkbook(fileName)
  sheets <- getSheets(wb)
  sheetNames <- toupper(names(sheets))
  
  truthFileIndex <- which(!is.na(match(sheetNames, "TRUTH")))
  if (truthFileIndex == 0) 
    stop("TRUTH table cannot be found in the dataset.")
  # truthTable <- readColumns(sheets[[truthFileIndex]], startColumn=1, endColumn = 3, startRow=1, endRow=NULL)
  truthTable <- read.xlsx2(fileName, truthFileIndex, colIndex = 1:3, colClasses = rep("numeric", 3))
  
  for (i in 1:3){
    truthTable[grep("^\\s*$", truthTable[ , i]), i] <- NA
  }
  
  naRows <- colSums(is.na(truthTable))
  if (max(naRows) > 0) {
    if (max(naRows) == min(naRows)) {
      truthTable <- truthTable[1:(nrow(truthTable) - max(naRows)), ]
    }
  }
  
  for (i in 1:2) {
    if (any((as.numeric(as.character(truthTable[, i]))) %% 1 != 0 )) {
      naLines <- which(!is.integer(as.numeric(as.character(truthTable[, i])))) + 1
      errorMsg <- paste0("There are non-integer values(s) for CaseID or LesionID at the line(s) ", paste(naLines, collapse = ", "), " in the TRUTH table.")
      stop(errorMsg)
    }
  }
  
  if (any(is.na(as.numeric(as.character(truthTable[, 3]))))) {
    naLines <- which(is.na(as.numeric(as.character(truthTable[, 3])))) + 1
    errorMsg <- paste0("There are non-numeric values(s) for weights at the line(s) ", paste(naLines, collapse = ", "), " in the TRUTH table.")
    stop(errorMsg)
  }
  
  caseID <- as.integer(truthTable[[1]])  # all 3 have same lenghts
  lesionID <- as.integer(truthTable[[2]])
  weights <- truthTable[[3]]
  
  normalCases <- sort(unique(caseID[lesionID == 0]))
  abnormalCases <- sort(unique(caseID[lesionID > 0]))
  allCases <- c(normalCases, abnormalCases)
  K1 <- length(normalCases)
  K2 <- length(abnormalCases)
  K <- (K1 + K2)
  
  if (anyDuplicated(cbind(caseID, lesionID))) {
    naLines <- which(duplicated(cbind(caseID, lesionID))) + 1
    errorMsg <- paste0("Line(s) ", paste(naLines, collapse = ", "), " in the TRUTH table are duplicated with previous line(s) .")
    stop(errorMsg)
  }
  
  nlFileIndex <- which(!is.na(match(sheetNames, c("FP", "NL"))))
  if (nlFileIndex == 0) 
    stop("FP table cannot be found in the dataset.")
  # NLTable <- readColumns(sheets[[nlFileIndex]], startColumn=1, endColumn = 4, startRow=1, endRow=NULL)
  NLTable <- read.xlsx2(fileName, nlFileIndex, colIndex = 1:4, colClasses = c(rep("character", 2), rep("numeric", 2)))
  
  for (i in 1:4){
    NLTable[grep("^\\s*$", NLTable[ , i]), i] <- NA
  }
  
  naRows <- colSums(is.na(NLTable))
  if (max(naRows) > 0) {
    if (max(naRows) == min(naRows)) {
      NLTable <- NLTable[1:(nrow(NLTable) - max(naRows)), ]
    }
  }
  
  for (i in 3:4) {
    if (any(is.na(as.numeric(as.character(NLTable[, i]))))) {
      naLines <- which(is.na(as.numeric(as.character(NLTable[, i])))) + 1
      errorMsg <- paste0("There are unavailable cell(s) at the line(s) ", paste(naLines, collapse = ", "), " in the FP table.")
      stop(errorMsg)
    }
  }
  
  NLReaderID <- as.character(NLTable[[1]])
  
  NLModalityID <- as.character(NLTable[[2]])
  
  NLCaseID <- NLTable[[3]]
  if (any(!(NLCaseID %in% caseID))) {
    naCases <- NLCaseID[which(!(NLCaseID %in% caseID))]
    errorMsg <- paste0("Case(s) ", paste(unique(naCases), collapse = ", "), " in the FP table cannot be found in TRUTH table.")
    stop(errorMsg)
  }
  NLRating <- NLTable[[4]]
  
  llFileIndex <- which(!is.na(match(sheetNames, c("TP", "LL"))))
  if (llFileIndex == 0) 
    stop("TP table cannot be found in the dataset.")
  # LLTable <- readColumns(sheets[[llFileIndex]], startColumn=1, endColumn = 5, startRow=1, endRow=NULL)
  LLTable <- read.xlsx2(fileName, llFileIndex, colIndex = 1:5, colClasses = c(rep("character", 2), rep("numeric", 3)))
  
  for (i in 1:5){
    LLTable[grep("^\\s*$", LLTable[ , i]), i] <- NA
  }
  
  naRows <- colSums(is.na(LLTable))
  if (max(naRows) > 0) {
    if (max(naRows) == min(naRows)) {
      LLTable <- LLTable[1:(nrow(LLTable) - max(naRows)), ]
    }
  }
  
  for (i in 3:5) {
    if (any(is.na(as.numeric(as.character(LLTable[, i]))))) {
      naLines <- which(is.na(as.numeric(as.character(LLTable[, i])))) + 1
      errorMsg <- paste0("There are unavailable cell(s) at the line(s) ", paste(naLines, collapse = ", "), " in the TP table.")
      stop(errorMsg)
    }
  }
  
  LLReaderID <- as.character(LLTable[[1]])
  
  LLModalityID <- as.character(LLTable[[2]])
  
  LLCaseID <- LLTable[[3]]
  LLLesionID <- LLTable[[4]]
  for (i in 1:nrow(LLTable)) {
    lineNum <- which((caseID == LLCaseID[i]) & (lesionID == LLLesionID[i]))
    if (!length(lineNum)) {
      errorMsg <- paste0("Modality ", LLTable[i, 2], " Reader(s) ", LLTable[i, 1], " Case(s) ", LLTable[i, 3], " Lesion(s) ", LLTable[i, 4], " cannot be found in TRUTH table .")
      stop(errorMsg)
    }
  }
  
  LLRating <- LLTable[[5]]
  
  if (anyDuplicated(LLTable[, 1:4])) {
    naLines <- which(duplicated(LLTable[, 1:4]))
    errorMsg <- paste0("Modality ", paste(LLTable[naLines, 2], collapse = ", "), " Reader(s) ", paste(LLTable[naLines, 1], collapse = ", "), " Case(s) ", paste(LLTable[naLines, 3], collapse = ", "), " Lesion(s) ", 
                       paste(LLTable[naLines, 4], collapse = ", "), " have multiple ratings in TP table .")
    stop(errorMsg)
  }
  
  lesionNum <- as.vector(table(caseID[caseID %in% abnormalCases]))
  # for (k2 in 1:length(abnormalCases)) { lesionNum[k2] <- sum(caseID == abnormalCases[k2]) }
  
  lesionWeight <- array(dim = c(length(abnormalCases), max(lesionNum)))
  lesionIDTable <- array(dim = c(length(abnormalCases), max(lesionNum)))
  
  for (k2 in 1:length(abnormalCases)) {
    k <- which(caseID == abnormalCases[k2])
    lesionIDTable[k2, ] <- c(sort(lesionID[k]), rep(UNINITIALIZED, max(lesionNum) - length(k)))
    if (all(weights[k] == 0)) {
      lesionWeight[k2, 1:length(k)] <- 1/lesionNum[k2]
    } else {
      lesionWeight[k2, ] <- c(weights[k][order(lesionID[k])], rep(UNINITIALIZED, max(lesionNum) - length(k)))
    }
  }
  
  modalityID <- as.character(sort(unique(c(NLModalityID, LLModalityID))))
  I <- length(modalityID)
  
  readerID <- as.character(sort(unique(c(NLReaderID, LLReaderID))))
  J <- length(readerID)
  
  maxNL <- 0
  for (i in modalityID) {
    for (j in readerID) {
      k <- (NLModalityID == i) & (NLReaderID == j)
      if ((sum(k) == 0)) 
        next
      maxNL <- max(maxNL, max(table(NLCaseID[k])))
    }
  }
  
  NL <- array(dim = c(I, J, K, maxNL))
  for (i in 1:I) {
    for (j in 1:J) {
      k <- (NLModalityID == modalityID[i]) & (NLReaderID == readerID[j])
      if ((sum(k) == 0)) 
        next
      caseNLTable <- table(NLCaseID[k])
      IDs <- as.numeric(unlist(attr(caseNLTable, "dimnames")))
      for (k1 in 1:length(IDs)) {
        for (el in 1:caseNLTable[k1]) {
          NL[i, j, which(IDs[k1] == allCases), el] <- NLRating[k][which(NLCaseID[k] == IDs[k1])][el]
        }
      }
    }
  }
  
  LL <- array(dim = c(I, J, K2, max(lesionNum)))
  for (i in 1:I) {
    for (j in 1:J) {
      k <- (LLModalityID == modalityID[i]) & (LLReaderID == readerID[j])
      if ((sum(k) == 0)) 
        next
      caseLLTable <- table(LLCaseID[k])
      IDs <- as.numeric(unlist(attr(caseLLTable, "dimnames")))
      for (k1 in 1:length(IDs)) {
        for (el in 1:caseLLTable[k1]) {
          LL[i, j, which(IDs[k1] == abnormalCases), which(LLLesionID[k][which(LLCaseID[k] == IDs[k1])][el] == lesionIDTable[which(IDs[k1] == abnormalCases), ])] <- LLRating[k][which(LLCaseID[k] == IDs[k1])][el]
        }
      }
    }
  }
  
  lesionWeight[is.na(lesionWeight)] <- UNINITIALIZED
  lesionIDTable[is.na(lesionIDTable)] <- UNINITIALIZED
  NL[is.na(NL)] <- UNINITIALIZED
  LL[is.na(LL)] <- UNINITIALIZED
  
  isROI <- TRUE
  for (i in 1:I) {
    for (j in 1:J) {
      if (any(NL[i, j, 1:K1, ] == UNINITIALIZED)) {
        isROI <- FALSE
        break
      }
      temp <- LL[i, j, , ] != UNINITIALIZED
      dim(temp) <- c(K2, max(lesionNum))
      if (!all(lesionNum == rowSums(temp))) {
        isROI <- FALSE
        break
      }
      temp <- NL[i, j, (K1 + 1):K, ] == UNINITIALIZED
      dim(temp) <- c(K2, maxNL)
      if (!all(lesionNum == rowSums(temp))) {
        isROI <- FALSE
        break
      }
    }
  }
  
  if ((max(table(caseID)) == 1) && (maxNL == 1) && (all((NL[, , (K1 + 1):K, ] == UNINITIALIZED))) && (all((NL[, , 1:K1, ] != UNINITIALIZED)))) {
    fileType <- "ROC"
  } else {
    if (isROI) {
      fileType <- "ROI"
    } else {
      fileType <- "FROC"
    }
  }
  
  return(list(NL = NL, LL = LL, lesionNum = lesionNum, lesionID = lesionIDTable, lesionWeight = lesionWeight, dataType = fileType, modalityID = modalityID, readerID = readerID))
} 
