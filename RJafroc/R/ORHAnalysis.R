#' Obuchowski-Rockette analysis with Hillis improvements
#' 
#' Performs Obuchowski-Rockette analysis with Hillis improvements for the specified dataset.
#' 
#' @param dataset See \link{DBMHAnalysis}.
#' @param fom See \link{DBMHAnalysis}.
#' @param alpha See \link{DBMHAnalysis}.
#' @param covEstMethod The method used to estimate the covariance matrix; can be \code{"Jackknife"}, \code{"Bootstrap"} or \code{"DeLong"}, 
#' the last assumes the Wilcoxon FOM is chosen, otherwise an error will result.
#' @param nBoots The number of bootstraps (default is 200), used only if the \code{"Bootstrap"} method is used to estimate 
#' the covariance matrix. 
#' @param option See \link{DBMHAnalysis}.
#' 
#' @return The return value is a list with following elements:
#' @return \item{fomArray}{Figures of merit array. See the return of \link{FigureOfMerit}.}
#' @return \item{msT}{Mean square of the figure of merit corresponding to the the treatment effect.}
#' @return \item{msTR}{Mean square of the figure of merit corresponding to the the treatment-reader effect.}
#' @return \item{varComp}{Obuchowski-Rockette variance component and covariance estimates.}
#' @return \item{fRRRC}{See \link{DBMHAnalysis}.}
#' @return \item{ddfRRRC}{See \link{DBMHAnalysis}.}
#' @return \item{pRRRC}{See \link{DBMHAnalysis}.}
#' @return \item{ciDiffTrtRRRC}{See \link{DBMHAnalysis}.}
#' @return \item{ciAvgRdrEachTrtRRRC}{See \link{DBMHAnalysis}.}
#' @return \item{fFRRC}{See \link{DBMHAnalysis}.}
#' @return \item{ddfFRRC}{See \link{DBMHAnalysis}.}
#' @return \item{pFRRC}{See \link{DBMHAnalysis}.}
#' @return \item{ciDiffTrtFRRC}{See \link{DBMHAnalysis}.}
#' @return \item{ciAvgRdrEachTrtFRRC}{See \link{DBMHAnalysis}.}
#' @return \item{ciDiffTrtEachRdr}{See \link{DBMHAnalysis}.}
#' @return \item{varCovEachRdr}{Obuchowski-Rockette Variance and Cov1 estimates for each reader.}
#' @return \item{fRRFC}{See \link{DBMHAnalysis}.}
#' @return \item{ddfRRFC}{See \link{DBMHAnalysis}.}
#' @return \item{pRRFC}{See \link{DBMHAnalysis}.}
#' @return \item{ciDiffTrtRRFC}{See \link{DBMHAnalysis}.}
#' @return \item{ciAvgRdrEachTrtRRFC}{See \link{DBMHAnalysis}.}
#' 
#' @examples
#' retOR <- ORHAnalysis(rocData, fom = "Wilcoxon")
#' print(retOR)
#' 
#' @export
#' @importFrom stats pf pt qt
#' 
ORHAnalysis <- function(dataset, fom = "wJAFROC", alpha = 0.05, covEstMethod = "Jackknife", nBoots = 200, option = "ALL") {
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
  
  if (!option %in% c("RRRC", "FRRC", "RRFC", "ALL")){
    errMsg <- sprintf("%s is not an available option.", option)
    stop(errMsg)
  }    
  
  if (I < 2) {
    stop("The analysis requires at least 2 modalities")
  }
  
  if (!covEstMethod %in% c("Jackknife", "Bootstrap", "DeLong")) {
    errMsg <- paste0(covEstMethod, " is not an allowed covariance estimation method.")
    stop(errMsg)
  }
  
  fomArray <- FigureOfMerit(dataset, fom)
  trMeans <- rowMeans(fomArray)
  fomMean <- mean(fomArray)
  
  ret <- EstimateVarCov(fomArray, NL, LL, lesionNum, lesionID, lesionWeight, maxNL, fom, covEstMethod, nBoots)
  var <- ret$var
  cov1 <- ret$cov1
  cov2 <- ret$cov2
  cov3 <- ret$cov3
  
  msT <- 0
  for (i in 1:I) {
    msT <- msT + (mean(fomArray[i, ]) - fomMean)^2
  }
  msT <- J * msT/(I - 1)
  
  msR <- 0
  for (j in 1:J) {
    msR <- msR + (mean(fomArray[, j]) - fomMean)^2
  }
  msR <- I * msR/(J - 1)
  
  msTR <- 0
  for (i in 1:I) {
    for (j in 1:J) {
      msTR <- msTR + (fomArray[i, j] - mean(fomArray[i, ]) - mean(fomArray[, j]) + fomMean)^2
    }
  }
  msTR <- msTR/((J - 1) * (I - 1))
  
  varTR <- msTR - var + cov1 + max(cov2 - cov3, 0)
  varR <- (msR - var - (I - 1) * cov1 + cov2 + (I - 1) * cov3 - varTR)/I
  varCovArray <- c(varR, varTR, cov1, cov2, cov3, var)
  nameArray <- c("Var(R)", "Var(T*R)", "COV1", "COV2", "COV3", "Var(Error)")
  varComp <- data.frame(varCov = varCovArray, row.names = nameArray)
  
  varSingle <- vector(length = I)
  cov2Single <- vector(length = I)
  for (i in 1:I) {
    fomSingle <- fomArray[i, ]
    nl <- NL[i, , , ]
    ll <- LL[i, , , ]
    dim(fomSingle) <- c(1, J)
    dim(nl) <- c(1, J, K, maxNL)
    dim(ll) <- c(1, J, K2, max(lesionNum))
    ret <- EstimateVarCov(fomSingle, nl, ll, lesionNum, lesionID, lesionWeight, maxNL, fom, covEstMethod, nBoots)
    varSingle[i] <- ret$var
    if (J > 1) {
      cov2Single[i] <- ret$cov2
    } else {
      cov2Single[i] <- 0
    }
  }
  
  varEchRder <- vector(length = J)
  cov1EchRder <- vector(length = J)
  for (j in 1:J) {
    fomSingle <- fomArray[, j]
    nl <- NL[, j, , ]
    ll <- LL[, j, , ]
    dim(fomSingle) <- c(I, 1)
    dim(nl) <- c(I, 1, K, maxNL)
    dim(ll) <- c(I, 1, K2, max(lesionNum))
    ret <- EstimateVarCov(fomSingle, nl, ll, lesionNum, lesionID, lesionWeight, maxNL, fom, covEstMethod, nBoots)
    varEchRder[j] <- ret$var
    cov1EchRder[j] <- ret$cov1
  }

  msRSingle <- array(0, dim = c(I))
  for (i in 1:I) {
    msRSingle[i] <- sum((fomArray[i, ] - trMeans[i])^2)/(J - 1)
  }
  
  diffTRMeans <- array(dim = choose(I, 2))
  diffTRName <- array(dim = choose(I, 2))
  ii <- 1
  for (i in 1:I) {
    if (i == I) 
      break
    for (ip in (i + 1):I) {
      diffTRMeans[ii] <- trMeans[i] - trMeans[ip]
      diffTRName[ii] <- paste(modalityID[i], modalityID[ip], sep = " - ")
      ii <- ii + 1
    }
  }
  
  msNum <- msT
  
  # ************ RRRC ****************
  if (option %in% c("RRRC", "ALL")) {
    if (J > 1) {
      msDenRRRC <- msTR + max(J * (cov2 - cov3), 0)
      fRRRC <- msNum/msDenRRRC
      ddfRRRC <- msDenRRRC^2/(msTR^2/((I - 1) * (J - 1)))
      pRRRC <- 1 - pf(fRRRC, I - 1, ddfRRRC)
      stdErrRRRC <- sqrt(2 * msDenRRRC/J)
      tStat <- vector()
      tPr <- vector()
      CIRRRC <- array(dim = c(length(diffTRMeans), 2))
      for (i in 1:length(diffTRMeans)) {
        tStat[i] <- diffTRMeans[i]/stdErrRRRC
        tPr[i] <- 2 * pt(abs(tStat[i]), ddfRRRC, lower.tail = FALSE) # critical correction, noted by user Lucy D'Agostino McGowan
        CIRRRC[i, ] <- sort(c(diffTRMeans[i] - qt(alpha/2, ddfRRRC) * stdErrRRRC, diffTRMeans[i] + qt(alpha/2, ddfRRRC) * stdErrRRRC))
      }
      ciDiffTrtRRRC <- data.frame(Treatment = diffTRName, Estimate = diffTRMeans, StdErr = rep(stdErrRRRC, choose(I, 2)), DF = rep(ddfRRRC, choose(I, 2)), t = tStat, p = tPr, CI = CIRRRC)
      colnames(ciDiffTrtRRRC) <- c("Treatment", "Estimate", "StdErr", "DF", "t", "Pr > t", "CI Lower", "CI Upper")
      
      dfSingleRRRC <- array(dim = I)
      msDenSingleRRRC <- array(dim = I)
      stdErrSingleRRRC <- array(dim = I)
      CISingleRRRC <- array(dim = c(I, 2))
      for (i in 1:I) {
        msDenSingleRRRC[i] <- msRSingle[i] + max(J * cov2Single[i], 0)
        dfSingleRRRC[i] <- msDenSingleRRRC[i]^2/msRSingle[i]^2 * (J - 1)
        stdErrSingleRRRC[i] <- sqrt(msDenSingleRRRC[i]/J)
        CISingleRRRC[i, ] <- sort(c(trMeans[i] - qt(alpha/2, dfSingleRRRC[i]) * stdErrSingleRRRC[i], trMeans[i] + qt(alpha/2, dfSingleRRRC[i]) * stdErrSingleRRRC[i]))
      }
      ciAvgRdrEachTrtRRRC <- data.frame(Treatment = modalityID, Area = trMeans, StdErr = stdErrSingleRRRC, DF = dfSingleRRRC, CI = CISingleRRRC, row.names = NULL)
      colnames(ciAvgRdrEachTrtRRRC) <- c("Treatment", "Area", "StdErr", "DF", "CI Lower", "CI Upper")
    } else {
      fRRRC <- NA
      ddfRRRC <- NA
      pRRRC <- NA
      ciDiffTrtRRRC <- NA
      ciAvgRdrEachTrtRRRC <- NA
    }
    if (option == "RRRC"){
      return(list(fomArray = fomArray, msT = msT, msTR = msTR, varComp = varComp, 
                  fRRRC = fRRRC, ddfRRRC = ddfRRRC, pRRRC = pRRRC, ciDiffTrtRRRC = ciDiffTrtRRRC, ciAvgRdrEachTrtRRRC = ciAvgRdrEachTrtRRRC))
    }
  }
  
  # ************ FRRC ****************
  if (option %in% c("FRRC", "ALL")) {
    if (J > 1) {
      msDenFRRC <- var - cov1 + (J - 1) * (cov2 - cov3)
    } else {
      msDenFRRC <- var - cov1
    }
    fFRRC <- msNum/msDenFRRC
    ddfFRRC <- Inf
    pFRRC <- 1 - pf(fFRRC, I - 1, ddfFRRC)
    stdErrFRRC <- sqrt(2 * msDenFRRC/J)
    tStat <- vector()
    tPr <- vector()
    CIFRRC <- array(dim = c(length(diffTRMeans), 2))
    for (i in 1:length(diffTRMeans)) {
      tStat[i] <- diffTRMeans[i]/stdErrFRRC
      tPr[i] <- 2 * pt(abs(tStat[i]), ddfFRRC, lower.tail = FALSE)  # critical correction, noted by user Lucy D'Agostino McGowan
      CIFRRC[i, ] <- sort(c(diffTRMeans[i] - qt(alpha/2, ddfFRRC) * stdErrFRRC, diffTRMeans[i] + qt(alpha/2, ddfFRRC) * stdErrFRRC))
    }
    ciDiffTrtFRRC <- data.frame(Treatment = diffTRName, Estimate = diffTRMeans, StdErr = rep(stdErrFRRC, choose(I, 2)), DF = rep(ddfFRRC, choose(I, 2)), t = tStat, p = tPr, CI = CIFRRC)
    colnames(ciDiffTrtFRRC) <- c("Treatment", "Estimate", "StdErr", "DF", "t", "Pr > t", "CI Lower", "CI Upper")
    
    dfSingleFRRC <- array(dim = I)
    msDenSingleFRRC <- array(dim = I)
    stdErrSingleFRRC <- array(dim = I)
    CISingleFRRC <- array(dim = c(I, 2))
    for (i in 1:I) {
      msDenSingleFRRC[i] <- varSingle[i] + (J - 1) * cov2Single[i]
      dfSingleFRRC[i] <- Inf
      stdErrSingleFRRC[i] <- sqrt(msDenSingleFRRC[i]/J)
      CISingleFRRC[i, ] <- sort(c(trMeans[i] - qt(alpha/2, dfSingleFRRC[i]) * stdErrSingleFRRC[i], trMeans[i] + qt(alpha/2, dfSingleFRRC[i]) * stdErrSingleFRRC[i]))
    }
    ciAvgRdrEachTrtFRRC <- data.frame(Treatment = modalityID, Area = trMeans, StdErr = stdErrSingleFRRC, DF = dfSingleFRRC, CI = CISingleFRRC, row.names = NULL)
    colnames(ciAvgRdrEachTrtFRRC) <- c("Treatment", "Area", "StdErr", "DF", "CI Lower", "CI Upper")
    
    diffTRMeansFRRC <- array(dim = c(J, choose(I, 2)))
    for (j in 1:J) {
      ii <- 1
      for (i in 1:I) {
        if (i == I) 
          break
        for (ip in (i + 1):I) {
          diffTRMeansFRRC[j, ii] <- fomArray[i, j] - fomArray[ip, j]
          ii <- ii + 1
        }
      }
    }
    
    diffTRMeansFRRC <- as.vector(t(diffTRMeansFRRC))
    stdErrFRRC <- sqrt(2 * (varEchRder - cov1EchRder))
    stdErrFRRC <- rep(stdErrFRRC, choose(I, 2))
    dim(stdErrFRRC) <- c(J, choose(I, 2))
    stdErrFRRC <- as.vector(t(stdErrFRRC))
    readerNames <- rep(readerID, choose(I, 2))
    dim(readerNames) <- c(J, choose(I, 2))
    readerNames <- as.vector(t(readerNames))
    trNames <- rep(diffTRName, J)
    dfReaderFRRC <- rep(Inf, length(stdErrFRRC))
    CIReaderFRRC <- array(dim = c(length(stdErrFRRC), 2))
    tStat <- vector()
    tPr <- vector()
    for (n in 1:length(stdErrFRRC)) {
      tStat[n] <- diffTRMeansFRRC[n]/stdErrFRRC[n]
      tPr[n] <- 2 * pt(abs(tStat[n]), dfReaderFRRC[n], lower.tail = FALSE)  # critical correction, noted by user Lucy D'Agostino McGowan
      CIReaderFRRC[n, ] <- sort(c(diffTRMeansFRRC[n] - qt(alpha/2, dfReaderFRRC[n]) * stdErrFRRC[n], diffTRMeansFRRC[n] + qt(alpha/2, dfReaderFRRC[n]) * stdErrFRRC[n]))
    }
    ciDiffTrtEachRdr <- data.frame(Reader = readerNames, Treatment = trNames, Estimate = diffTRMeansFRRC, StdErr = stdErrFRRC, DF = dfReaderFRRC, t = tStat, p = tPr, CI = CIReaderFRRC)
    colnames(ciDiffTrtEachRdr) <- c("Reader", "Treatment", "Estimate", "StdErr", "DF", "t", "Pr > t", "CI Lower", "CI Upper")
    
    varCovEachRdr <- data.frame(readerID, varEchRder, cov1EchRder)
    colnames(varCovEachRdr) <- c("Reader", "Var", "Cov1")
    if (option == "FRRC"){
      return(list(fomArray = fomArray, msT = msT, msTR = msTR, varComp = varComp, 
                  fFRRC = fFRRC, ddfFRRC = ddfFRRC, pFRRC = pFRRC, ciDiffTrtFRRC = ciDiffTrtFRRC, ciAvgRdrEachTrtFRRC = ciAvgRdrEachTrtFRRC, ciDiffTrtEachRdr = ciDiffTrtEachRdr, varCovEachRdr = varCovEachRdr
      ))
    }
  }
  
  # ************ RRFC ****************
  if (option %in% c("RRFC", "ALL")) {
    if (J > 1) {
      msDenRRFC <- msTR
      fRRFC <- msNum/msDenRRFC
      ddfRRFC <- ((I - 1) * (J - 1))
      pRRFC <- 1 - pf(fRRFC, I - 1, ddfRRFC)
      stdErrRRFC <- sqrt(2 * msDenRRFC/J)
      tStat <- vector()
      tPr <- vector()
      CIRRFC <- array(dim = c(length(diffTRMeans), 2))
      for (i in 1:length(diffTRMeans)) {
        tStat[i] <- diffTRMeans[i]/stdErrRRFC
        tPr[i] <- 2 * pt(abs(tStat[i]), ddfRRFC, lower.tail = FALSE)  # critical correction, noted by user Lucy D'Agostino McGowan
        CIRRFC[i, ] <- sort(c(diffTRMeans[i] - qt(alpha/2, ddfRRFC) * stdErrRRFC, diffTRMeans[i] + qt(alpha/2, ddfRRFC) * stdErrRRFC))
      }
      ciDiffTrtRRFC <- data.frame(Treatment = diffTRName, Estimate = diffTRMeans, StdErr = rep(stdErrRRFC, choose(I, 2)), DF = rep(ddfRRFC, choose(I, 2)), t = tStat, p = tPr, CI = CIRRFC)
      colnames(ciDiffTrtRRFC) <- c("Treatment", "Estimate", "StdErr", "DF", "t", "Pr > t", "CI Lower", "CI Upper")
      
      dfSingleRRFC <- array(dim = I)
      msDenSingleRRFC <- array(dim = I)
      stdErrSingleRRFC <- array(dim = I)
      CISingleRRFC <- array(dim = c(I, 2))
      for (i in 1:I) {
        msDenSingleRRFC[i] <- msRSingle[i]
        dfSingleRRFC[i] <- (J - 1)
        stdErrSingleRRFC[i] <- sqrt(msDenSingleRRFC[i]/J)
        CISingleRRFC[i, ] <- sort(c(trMeans[i] - qt(alpha/2, dfSingleRRFC[i]) * stdErrSingleRRFC[i], trMeans[i] + qt(alpha/2, dfSingleRRFC[i]) * stdErrSingleRRFC[i]))
      }
      ciAvgRdrEachTrtRRFC <- data.frame(Treatment = modalityID, Area = trMeans, StdErr = stdErrSingleRRFC, DF = dfSingleRRFC, CI = CISingleRRFC, row.names = NULL)
      colnames(ciAvgRdrEachTrtRRFC) <- c("Treatment", "Area", "StdErr", "DF", "CI Lower", "CI Upper")
    } else {
      fRRFC <- NA
      ddfRRFC <- NA
      pRRFC <- NA
      ciDiffTrtRRFC <- NA
      ciAvgRdrEachTrtRRFC <- NA
    }
    if (option == "RRFC"){
      return(list(fomArray = fomArray, msT = msT, msTR = msTR, varComp = varComp, 
                  fRRFC = fRRFC, ddfRRFC = ddfRRFC, pRRFC = pRRFC, ciDiffTrtRRFC = ciDiffTrtRRFC, ciAvgRdrEachTrtRRFC = ciAvgRdrEachTrtRRFC))
    }
  }
  
  return(list(fomArray = fomArray, msT = msT, msTR = msTR, varComp = varComp, 
              fRRRC = fRRRC, ddfRRRC = ddfRRRC, pRRRC = pRRRC, ciDiffTrtRRRC = ciDiffTrtRRRC, ciAvgRdrEachTrtRRRC = ciAvgRdrEachTrtRRRC, 
              fFRRC = fFRRC, ddfFRRC = ddfFRRC, pFRRC = pFRRC, ciDiffTrtFRRC = ciDiffTrtFRRC, ciAvgRdrEachTrtFRRC = ciAvgRdrEachTrtFRRC, ciDiffTrtEachRdr = ciDiffTrtEachRdr, varCovEachRdr = varCovEachRdr, 
              fRRFC = fRRFC, ddfRRFC = ddfRRFC, pRRFC = pRRFC, ciDiffTrtRRFC = ciDiffTrtRRFC, ciAvgRdrEachTrtRRFC = ciAvgRdrEachTrtRRFC))
} 
