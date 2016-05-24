#' Convert FROC dataset
#' 
#' Convert an FROC dataset to a highest rating inferred ROC dataset
#' 
#' @param dataset The dataset to be converted, see \link{RJafroc-package}.
#' 
#' @return An ROC dataset, where each case is represented by the highest rating of the case in FROC dataset.
#'  
#'
#' @examples
#' FROC2HrROC(dataset = frocData)
#' 
#' @export
#'  

FROC2HrROC <- function (dataset){
  UNINITIALIZED <- RJafrocEnv$UNINITIALIZED
  NL <- dataset$NL
  LL <- dataset$LL
  I <- length(dataset$modalityID)
  J <- length(dataset$readerID)
  K <- dim(NL)[3]
  K2 <- dim(LL)[3]
  K1 <- K - K2 
  modalityID <- dataset$modalityID
  readerID <- dataset$readerID
  maxNL <- dim(NL)[4] 
  lesionNum <- dataset$lesionNum
  
  NL <- apply(NL, c(1, 2, 3), max)
  LL <- apply(LL, c(1, 2, 3), max)
  
  dim(NL) <- c(dim(NL), 1)
  dim(LL) <- c(dim(LL), 1)
  
  LLTmp <- array(dim = c(I, J, K2, 2))
  LLTmp[ , , , 1] <- NL[ , , (K1 + 1):K, ]
  LLTmp[ , , , 2] <- LL
  LL <- apply(LLTmp, c(1, 2, 3), max)
  dim(LL) <- c(dim(LL), 1)
  NL[NL == UNINITIALIZED] <- -2000
  LL[LL == UNINITIALIZED] <- -2000
  
  NL[ , , (K1 + 1):K, ] <- UNINITIALIZED
  
  lesionNum <- rep(1, times = length(lesionNum))
  lesionID <- lesionNum
  dim(lesionID) <- c(K2, 1)
  lesionWeight <- lesionID 
  dataset$NL <- NL
  dataset$LL <- LL
  dataset$lesionNum <- lesionNum
  dataset$lesionID <- lesionID
  dataset$lesionWeight <- lesionWeight
  dataset$dataType <- "ROC"
  
  return (dataset)
}