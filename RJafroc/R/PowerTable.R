#' Calculate power table, different combinations, of J and K for desired power for ROC studies.
#' 
#' Calculate required sample size for the specified dataset with given significance level, effect size and desired power
#' 
#' @param dataset The dataset to be analyzed, see \link{RJafroc-package}.
#' @param alpha The significantce level.
#' @param effectSize The effect size to be used in the calculation.
#' @param desiredPower The desired statistical power.
#' @param randomOption The random option. It can be \code{"ALL"}, \code{"READERS"} or \code{"CASES"}, which indicate predictions for (1) random readers and random cases, (2) random readers only and 
#' (3) random cases only.
#' 
#' @return The return is a data frame containing following three columns.
#' @return \item{numReaders}{The number of readers.}  
#' @return \item{numCases}{The number of cases.}
#' @return \item{power}{The statistical power for the number of readers and cases combination.}
#' 
#' @examples
#' retDbm <- DBMHAnalysis(data = rocData, fom = "Wilcoxon")                     
#' effectSize <- retDbm$ciDiffTrtRRRC$Estimate
#' PowerTable(dataset = rocData, effectSize = effectSize)
#' 
#' @export
#' 
PowerTable <- function(dataset, alpha = 0.05, effectSize = 0.05, desiredPower = 0.8, randomOption = "ALL") {
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
  
  dim(NL) <- c(I, J, K, maxNL)
  dim(LL) <- c(I, J, K2, max(lesionNum))
  
  if (J < 2) 
    stop("The dataset must have at least 2 readers.")
  if (dataType == "ROI") 
    stop("Cannot calculate sample size for ROI data.")
  fomArray <- array(dim = c(I, J))
  for (i in 1:I) {
    for (j in 1:J) {
      nl <- NL[i, j, , ]
      ll <- LL[i, j, , ]
      dim(nl) <- c(K, maxNL)
      dim(ll) <- c(K2, max(lesionNum))
      fomArray[i, j] <- MyFOM(nl, ll, lesionNum, lesionID, lesionWeight, maxNL, "Wilcoxon")
    }
  }
  
  fomMean <- mean(fomArray)
  
  msTR <- 0
  for (i in 1:I) {
    for (j in 1:J) {
      msTR <- msTR + (fomArray[i, j] - mean(fomArray[i, ]) - mean(fomArray[, j]) + fomMean)^2
    }
  }
  msTR <- msTR/((J - 1) * (I - 1))
  
  ret <- EstimateVarCov(fomArray, NL, LL, lesionNum, lesionID, lesionWeight, maxNL, "HrAuc", "Jackknife")
  var <- ret$var
  cov1 <- ret$cov1
  cov2 <- ret$cov2
  cov3 <- ret$cov3
  
  varYTR <- msTR - var + cov1 + max(cov2 - cov3, 0)
  varYTC <- max(cov2 - cov3, 0) * K
  varYRC <- max(cov1 - cov3, 0) * K
  varYEps <- (var - cov2) * K - varYRC
  
  nCases <- 2000
  j <- 1
  randomSampleSize <- NULL
  while (nCases >= 20) {
    j <- j + 1
    if (j > 10) 
      break
    ret <- SampleSizeGivenJ(j, varYTR, varYTC, varYEps, alpha = alpha, effectSize = effectSize, desiredPower = desiredPower, randomOption = randomOption)
    nCases <- ret$K
    power <- ret$power
    if (nCases > 2000) {
      randomSampleSize <- rbind(randomSampleSize, c(j, ">2000", NA))
    } else if (nCases < 20) {
      randomSampleSize <- rbind(randomSampleSize, c(j, "<20", NA))
    } else {
      randomSampleSize <- rbind(randomSampleSize, c(j, nCases, signif(power, 3)))
    }
  }
  randomSampleSize <- data.frame(numReaders = randomSampleSize[, 1], numCases = randomSampleSize[, 2], power = randomSampleSize[, 3])
  return(randomSampleSize)
} 
