#' @importFrom tools file_ext
#' @importFrom stringr str_trim
ReadLrc <- function(fileName) {
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
  fileExt <- file_ext(fileName)
  if (fileExt != "lrc") {
    stop("The extension of LRC file name must be *.lrc")
  }
  whiteSpaces <- c("", " ", "\t")
  
  fileLines <- readLines(fileName)
  fileLines <- str_trim(fileLines)
  lineNum <- 2
  
  J <- 1
  readerID <- fileLines[lineNum]
  lineNum <- 3
  modalityID <- unlist(strsplit(fileLines[lineNum], split = "\""))
  iRmv <- NULL
  ifRmv <- FALSE
  for (i in 1:length(modalityID)) {
    modalityChar <- unlist(strsplit(modalityID[i], split = ""))
    if (all(modalityChar %in% whiteSpaces)) {
      iRmv <- c(iRmv, i)
      ifRmv <- TRUE
    }
  }
  if (ifRmv) {
    modalityID <- modalityID[-iRmv]
  }
  I <- length(modalityID)
  
  lineNum <- 4
  largeOrSmall <- splitWhiteSpaces(fileLines[lineNum])
  if (length(largeOrSmall) != I) {
    stop("The numbers of elements in the third and fourth lines does not match.")
  }
  
  lineNum <- 5
  NL <- NA
  LL <- NA
  while (TRUE) {
    K1 <- NA
    K1Temp <- 0
    NLSnglR <- NULL
    while (fileLines[lineNum] != "*") {
      K1Temp <- K1Temp + 1
      ratings <- splitWhiteSpaces(fileLines[lineNum])
      if (length(ratings) < I) {
        errMsg <- sprintf("The number of ratings in line %d is less than the number of maodalities.", lineNum)
        stop(errMsg)
      }
      NLSnglR <- rbind(NLSnglR, ratings[1:I])
      lineNum <- lineNum + 1
    }
    if (is.na(K1)) {
      K1 <- K1Temp
    } else if (K1 != K1Temp) {
      errMsg <- sprintf("The number of non-diseased cases for reader %d does not match that for previous readers.", J)
      stop(errMsg)
    }
    NLTemp <- array(dim = c(I, J, K1, 1))
    if (is.null(dim(NL))) {
      NLTemp[, , , 1] <- as.numeric(t(NLSnglR))
      NL <- NLTemp
    } else {
      NLTemp[, 1:(J - 1), , ] <- NL
      NLTemp[, J, , 1] <- as.numeric(t(NLSnglR))
      NL <- NLTemp
    }
    
    K2 <- NA
    K2Temp <- 0
    LLSnglR <- NULL
    lineNum <- lineNum + 1
    while (fileLines[lineNum] != "*") {
      K2Temp <- K2Temp + 1
      ratings <- splitWhiteSpaces(fileLines[lineNum])
      if (length(ratings) < I) {
        errMsg <- sprintf("The number of ratings in line %d is less than the number of maodalities.", lineNum)
        stop(errMsg)
      }
      LLSnglR <- rbind(LLSnglR, ratings[1:I])
      lineNum <- lineNum + 1
    }
    if (is.na(K2)) {
      K2 <- K2Temp
    } else if (K2 != K2Temp) {
      errMsg <- sprintf("The number of diseased cases for reader %d does not match that for previous readers.", J)
      stop(errMsg)
    }
    LLTemp <- array(dim = c(I, J, K2, 1))
    if (is.null(dim(LL))) {
      LLTemp[, , , 1] <- as.numeric(t(LLSnglR))
      LL <- LLTemp
    } else {
      LLTemp[, 1:(J - 1), , ] <- LL
      LLTemp[, J, , 1] <- as.numeric(t(LLSnglR))
      LL <- LLTemp
    }
    lineNum <- lineNum + 1
    if (fileLines[lineNum] != "#") {
      J <- J + 1
      readerID <- c(readerID, fileLines[lineNum])
      lineNum <- lineNum + 1
    } else {
      break
    }
  }
  for (i in 1:I) {
    if (toupper(largeOrSmall[i]) %in% c("S", "SMALL")) {
      NL[i, , , ] <- max(NL[i, , , ]) + 1 - NL[i, , , ]
    } else if (!toupper(largeOrSmall[i]) %in% c("L", "LARGE")) {
      stop("The confidence indicater in line 4 must be \"L\", \"LARGE\", \"S\", or \"SMALL\".")
    }
  }
  K <- K1 + K2
  NLTemp <- array(UNINITIALIZED, dim = c(I, J, K, 1))
  NLTemp[, , 1:K1, ] <- NL
  NL <- NLTemp
  lesionNum <- rep(1, K2)
  lesionID <- array(1, dim = c(K2, 1))
  lesionWeight <- lesionID
  maxNL <- 1
  dataType <- "ROC"
  return(list(NL = NL, LL = LL, lesionNum = lesionNum, lesionID = lesionID, lesionWeight = lesionWeight, dataType = dataType, modalityID = modalityID, readerID = readerID))
}

splitWhiteSpaces <- function(string) {
  whiteSpaces <- c("", " ", "\t")
  string <- unlist(strsplit(string, split = " |\t"))
  string <- string[!string %in% whiteSpaces]
  return(string)
} 
