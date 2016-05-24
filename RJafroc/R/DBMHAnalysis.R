#' DBM analysis with Hillis improvements
#' 
#' Performs Dorfman-Berbaum-Metz significance testing, with Hillis improvements, for the specified dataset.
#' 
#' @param dataset The dataset to be analyzed, see \link{RJafroc-package}.
#' @param fom The figure of merit to be used in the analysis, default is \code{"wJAFROC"}, see \link{FigureOfMerit}.
#' @param alpha The significance level of the test of the null hypothesis that all treatment effects are zero (default alpha is 0.05).
#' @param option The analysis option: it can be \code{"RRRC"}, \code{"FRRC"}, \code{"RRFC"} or \code{"ALL"} (the default), corresponding 
#' random readers and random cases, fixed readers and random cases and random readers and fixed cases, respectively; this option 
#' indicates which factors are treated as random and/or fixed factors in the analysis.
#' 
#' @return The return value is a list with following elements:
#' @return \item{fomArray}{The figure of merit array of each reader and modality.}
#' @return \item{anovaY}{The ANOVA table of the pseudovalues.}
#' @return \item{anovaYi}{The ANOVA table of the pseudovalues for each modality.}
#' @return \item{varComp}{The DBM variance components estimates.}
#' @return \item{fRRRC}{The F statistic for testing the null hypothesis, for the RRRC condition.}
#' @return \item{ddfRRRC}{The denominator degrees of freedom of the F statistic, for the RRRC condition.}
#' @return \item{pRRRC}{The p-value of the significance test of the NH for the RRRC condition.}
#' @return \item{ciDiffTrtRRRC}{The confidence intervals and related tests for the FOM differences between pairs of modalities, 
#' for the RRRC condition.}
#' @return \item{ciAvgRdrEachTrtRRRC}{The confidence intervals and related tests for rdr. avg. FOM in each modality, for the RRRC condition.}
#' @return \item{fFRRC}{The F statistic for testing the null hypothesis, for the FRRC condition.}
#' @return \item{ddfFRRC}{The denominator degrees of freedom of the FRRC F statistic.}
#' @return \item{pFRRC}{The p-value of the significance test of the NH, for the FRRC condition.}
#' @return \item{ciDiffTrtFRRC}{The confidence intervals and related tests for the FOM differences between pairs of modalities, 
#' for the FRRC condition.}
#' @return \item{ciAvgRdrEachTrtFRRC}{The confidence intervals and related tests for rdr. avg. FOM in each modality, for the FRRC condition.}
#' @return \item{ssAnovaEachRdr}{The sum of squares table of the ANOVA of the pseudovalues for each reader (based on the data only 
#' for the specified reader).}
#' @return \item{msAnovaEachRdr}{The mean squares table of the ANOVA of the pseudovalues for each reader (based on the data only 
#' for the specified reader).}
#' @return \item{ciDiffTrtEachRdr}{The confidence intervals and related tests of the FOM differences between pairs of modalities for each reader.}
#' @return \item{fRRFC}{The F statistic for testing the null hypothesis, for the RRFC condition.}
#' @return \item{ddfRRFC}{The denominator degrees of freedom of the F statistic, for the RRFC condition.}
#' @return \item{pRRFC}{The p-value of the significance test of the NH, for the RRFC condition.}
#' @return \item{ciDiffTrtRRFC}{The confidence intervals and related tests for the FOM differences between pairs of modalities, 
#' for the RRFC condition.}
#' @return \item{ciAvgRdrEachTrtRRFC}{The confidence intervals and related tests for rdr. avg. FOM in each modality, for the RRFC condition.}
#'  
#' @examples
#' retDbmRoc  <- DBMHAnalysis(rocData, fom = "Wilcoxon") 
#' 
#' \dontrun{
#' retDbmwJAFROC  <- DBMHAnalysis(frocData) # default is weighted JAFROC
#' 
#' retDbmHrAuc  <- DBMHAnalysis(frocData, fom = "HrAuc") 
#' print(retDbmHrAuc) 
#' 
#' retDbmSongA1  <- DBMHAnalysis(frocData, fom = "SongA1") 
#' print(retDbmSongA1)
#' 
#' retDbmSongA2  <- DBMHAnalysis(frocData, fom = "SongA2") 
#' print(retDbmSongA2)
#' 
#' retDbmwJafroc1  <- DBMHAnalysis(frocData, fom = "wJAFROC1")
#' print(retDbmwJafroc1)
#'  
#' retDbmJafroc1  <- DBMHAnalysis(frocData, fom = "JAFROC1")
#' print(retDbmJafroc1)
#'  
#' retDbmJAFROC  <- DBMHAnalysis(frocData, fom = "JAFROC") 
#' print(retDbmJAFROC)
#'  
#' retDbmROI  <- DBMHAnalysis(roiData, fom = "ROI")
#' print(retDbmROI)
#' }
#' 
#' @export
#' @importFrom stats pf pt qt
DBMHAnalysis <- function(dataset, fom = "wJAFROC", alpha = 0.05, option = "ALL") {
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
  
  if (!option %in% c("RRRC", "FRRC", "RRFC", "ALL")){
    errMsg <- sprintf("%s is not an available option.", option)
    stop(errMsg)
  }    
  
  if (I < 2) {
    stop("The analysis requires at least 2 modalities; consult Dr. Chakraborty for single-modality analysis")
  }
  
  fomArray <- FigureOfMerit(dataset, fom)
  trMeans <- rowMeans(fomArray)
  
  if (fom %in% c("MaxNLF", "ExpTrnsfmSp", "HrSp")) {
    jkFOMArray <- array(dim = c(I, J, K1))
    pseudoValues <- array(dim = c(I, J, K1))
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K1) {
          nl <- NL[i, j, -k, ]
          ll <- LL[i, j, , ]
          dim(nl) <- c(K - 1, maxNL)
          dim(ll) <- c(K2, max(lesionNum))
          jkFOMArray[i, j, k] <- MyFOM(nl, ll, lesionNum, lesionID, lesionWeight, maxNL, fom)
          pseudoValues[i, j, k] <- fomArray[i, j] * K1 - jkFOMArray[i, j, k] * (K1 - 1)
        }
        pseudoValues[i, j, ] <- pseudoValues[i, j, ] + (fomArray[i, j] - mean(pseudoValues[i, j, ]))
      }
    }
  } else if (fom %in% c("MaxLLF", "HrSe")) {
    jkFOMArray <- array(dim = c(I, J, K2))
    pseudoValues <- array(dim = c(I, J, K2))
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K2) {
          nl <- NL[i, j, -(k + K1), ]
          ll <- LL[i, j, -k, ]
          dim(nl) <- c(K - 1, maxNL)
          dim(ll) <- c(K2 - 1, max(lesionNum))
          jkFOMArray[i, j, k] <- MyFOM(nl, ll, lesionNum[-k], lesionID[-k, ], lesionWeight[-k, ], maxNL, fom)
          pseudoValues[i, j, k] <- fomArray[i, j] * K2 - jkFOMArray[i, j, k] * (K2 - 1)
        }
        pseudoValues[i, j, ] <- pseudoValues[i, j, ] + (fomArray[i, j] - mean(pseudoValues[i, j, ]))
      }
    }
  } else {
    jkFOMArray <- array(dim = c(I, J, K))
    pseudoValues <- array(dim = c(I, J, K))
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K) {
          if (k <= K1) {
            nl <- NL[i, j, -k, ]
            ll <- LL[i, j, , ]
            dim(nl) <- c(K - 1, maxNL)
            dim(ll) <- c(K2, max(lesionNum))
            jkFOMArray[i, j, k] <- MyFOM(nl, ll, lesionNum, lesionID, lesionWeight, maxNL, fom)
          } else {
            nl <- NL[i, j, -k, ]
            ll <- LL[i, j, -(k - K1), ]
            dim(nl) <- c(K - 1, maxNL)
            dim(ll) <- c(K2 - 1, max(lesionNum))
            jkFOMArray[i, j, k] <- MyFOM(nl, ll, lesionNum[-(k - K1)], lesionID[-(k - K1), ], lesionWeight[-(k - K1), ], maxNL, fom)
          }
          pseudoValues[i, j, k] <- fomArray[i, j] * K - jkFOMArray[i, j, k] * (K - 1)
        }
        pseudoValues[i, j, ] <- pseudoValues[i, j, ] + (fomArray[i, j] - mean(pseudoValues[i, j, ]))
      }
    }
  }
  K <- length(pseudoValues[1, 1, ])
  
  msT <- 0
  for (i in 1:I) {
    msT <- msT + (mean(pseudoValues[i, , ]) - mean(pseudoValues))^2
  }
  msT <- msT * K * J/(I - 1)
  
  
  msR <- 0
  for (j in 1:J) {
    msR <- msR + (mean(pseudoValues[, j, ]) - mean(pseudoValues))^2
  }
  msR <- msR * K * I/(J - 1)
  
  
  msC <- 0
  for (k in 1:K) {
    msC <- msC + (mean(pseudoValues[, , k]) - mean(pseudoValues))^2
  }
  msC <- msC * I * J/(K - 1)
  
  
  msTR <- 0
  for (i in 1:I) {
    for (j in 1:J) {
      msTR <- msTR + (mean(pseudoValues[i, j, ]) - mean(pseudoValues[i, , ]) - mean(pseudoValues[, j, ]) + mean(pseudoValues))^2
    }
  }
  msTR <- msTR * K/((I - 1) * (J - 1))
  
  
  msTC <- 0
  for (i in 1:I) {
    for (k in 1:K) {
      msTC <- msTC + (mean(pseudoValues[i, , k]) - mean(pseudoValues[i, , ]) - mean(pseudoValues[, , k]) + mean(pseudoValues))^2
    }
  }
  msTC <- msTC * J/((I - 1) * (K - 1))
  
  
  msRC <- 0
  for (j in 1:J) {
    for (k in 1:K) {
      msRC <- msRC + (mean(pseudoValues[, j, k]) - mean(pseudoValues[, j, ]) - mean(pseudoValues[, , k]) + mean(pseudoValues))^2
    }
  }
  msRC <- msRC * I/((J - 1) * (K - 1))
  
  msTRC <- 0
  for (i in 1:I) {
    for (j in 1:J) {
      for (k in 1:K) {
        msTRC <- msTRC + (pseudoValues[i, j, k] - mean(pseudoValues[i, j, ]) - mean(pseudoValues[i, , k]) - mean(pseudoValues[, j, k]) + 
                            mean(pseudoValues[i, , ]) + mean(pseudoValues[, j, ]) + mean(pseudoValues[, , k]) - mean(pseudoValues))^2
      }
    }
  }
  msTRC <- msTRC/((I - 1) * (J - 1) * (K - 1))
  
  varR <- (msR - msTR - msRC + msTRC)/(I * K)
  varC <- (msC - msTC - msRC + msTRC)/(I * J)
  varTR <- (msTR - msTRC)/K
  varTC <- (msTC - msTRC)/J
  varRC <- (msRC - msTRC)/I
  varErr <- msTRC
  varComp <- c(varR, varC, varTR, varTC, varRC, varErr)
  varCompName <- c("Var(R)", "Var(C)", "Var(T*R)", "Var(T*C)", "Var(R*C)", "Var(Error)")
  varComp <- data.frame(varComp, row.names = varCompName)
  
  msArray <- c(msT, msR, msC, msTR, msTC, msRC, msTRC)
  dfArray <- c(I - 1, J - 1, K - 1, (I - 1) * (J - 1), (I - 1) * (K - 1), (J - 1) * (K - 1), (I - 1) * (J - 1) * (K - 1))
  ssArray <- msArray * dfArray
  msArray <- c(msArray, NA)
  dfArray <- c(dfArray, sum(dfArray))
  ssArray <- c(ssArray, sum(ssArray))
  sourceArray <- c("T", "R", "C", "TR", "TC", "RC", "TRC", "Total")
  anovaY <- data.frame(Source = sourceArray, SS = ssArray, DF = dfArray, MS = msArray)
  
  msRSingle <- array(0, dim = c(I))
  msCSingle <- array(0, dim = c(I))
  msRCSingle <- array(0, dim = c(I))
  for (i in 1:I) {
    for (j in 1:J) {
      msRSingle[i] <- msRSingle[i] + (mean(pseudoValues[i, j, ]) - mean(pseudoValues[i, , ]))^2
    }
    msRSingle[i] <- msRSingle[i] * K/(J - 1)
    
    
    for (k in 1:K) {
      msCSingle[i] <- msCSingle[i] + (mean(pseudoValues[i, , k]) - mean(pseudoValues[i, , ]))^2
    }
    msCSingle[i] <- msCSingle[i] * J/(K - 1)
    
    for (j in 1:J) {
      for (k in 1:K) {
        msRCSingle[i] <- msRCSingle[i] + (mean(pseudoValues[i, j, k]) - mean(pseudoValues[i, j, ]) - mean(pseudoValues[i, , k]) + mean(pseudoValues[i, , ]))^2
      }
    }
    msRCSingle[i] <- msRCSingle[i]/((J - 1) * (K - 1))
  }
  sourceArraySingle <- c("R", "C", "RC")
  dfArraySingle <- c(J - 1, K - 1, (J - 1) * (K - 1))
  msArraySingle <- t(cbind(msRSingle, msCSingle, msRCSingle))
  msSingleTable <- data.frame(sourceArraySingle, dfArraySingle, msArraySingle, row.names = NULL)
  colnames(msSingleTable) <- c("Source", "DF", modalityID)
  
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
  
  if (option %in% c("RRRC", "ALL")) {
    # ************ RRRC ****************
    if (J > 1) {
      msDenRRRC <- msTR + max(msTC - msTRC, 0)
      fRRRC <- msNum/msDenRRRC
      ddfRRRC <- msDenRRRC^2/(msTR^2/((I - 1) * (J - 1)))
      pRRRC <- 1 - pf(fRRRC, I - 1, ddfRRRC)
      stdErrRRRC <- sqrt(2 * msDenRRRC/J/K)
      tStat <- vector()
      tPr <- vector()
      CIRRRC <- array(dim = c(length(diffTRMeans), 2))
      for (i in 1:length(diffTRMeans)) {
        tStat[i] <- diffTRMeans[i]/stdErrRRRC
        tPr[i] <- 2 * pt(abs(tStat[i]), ddfRRRC, lower.tail = FALSE)  # critical correction, noted by user Lucy D'Agostino McGowan
        CIRRRC[i, ] <- sort(c(diffTRMeans[i] - qt(alpha/2, ddfRRRC) * stdErrRRRC, diffTRMeans[i] + qt(alpha/2, ddfRRRC) * stdErrRRRC))
        
      }
      ciDiffTrtRRRC <- data.frame(Treatment = diffTRName, Estimate = diffTRMeans, StdErr = rep(stdErrRRRC, choose(I, 2)), DF = rep(ddfRRRC, choose(I, 2)), t = tStat, p = tPr, CI = CIRRRC)
      colnames(ciDiffTrtRRRC) <- c("Treatment", "Estimate", "StdErr", "DF", "t", "Pr > t", "CI Lower", "CI Upper")
      
      dfSingleRRRC <- array(dim = I)
      msDenSingleRRRC <- array(dim = I)
      stdErrSingleRRRC <- array(dim = I)
      CISingleRRRC <- array(dim = c(I, 2))
      for (i in 1:I) {
        msDenSingleRRRC[i] <- msRSingle[i] + max(msCSingle[i] - msRCSingle[i], 0)
        dfSingleRRRC[i] <- msDenSingleRRRC[i]^2/msRSingle[i]^2 * (J - 1)
        stdErrSingleRRRC[i] <- sqrt(msDenSingleRRRC[i]/J/K)
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
    if (option == "RRRC")
      return(list(fomArray = fomArray, anovaY = anovaY, anovaYi = msSingleTable, varComp = varComp, 
                  fRRRC = fRRRC, ddfRRRC = ddfRRRC, pRRRC = pRRRC, ciDiffTrtRRRC = ciDiffTrtRRRC, ciAvgRdrEachTrtRRRC = ciAvgRdrEachTrtRRRC))
  }
  
  if (option %in% c("FRRC", "ALL")) {
    # ************ FRRC ****************
    msDenFRRC <- msTC
    fFRRC <- msNum/msDenFRRC
    ddfFRRC <- (I - 1) * (K - 1)
    pFRRC <- 1 - pf(fFRRC, I - 1, ddfFRRC)
    stdErrFRRC <- sqrt(2 * msDenFRRC/J/K)
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
      msDenSingleFRRC[i] <- msCSingle[i]
      dfSingleFRRC[i] <- (K - 1)
      stdErrSingleFRRC[i] <- sqrt(msDenSingleFRRC[i]/J/K)
      CISingleFRRC[i, ] <- sort(c(trMeans[i] - qt(alpha/2, dfSingleFRRC[i]) * stdErrSingleFRRC[i], trMeans[i] + qt(alpha/2, dfSingleFRRC[i]) * stdErrSingleFRRC[i]))
    }
    ciAvgRdrEachTrtFRRC <- data.frame(Treatment = modalityID, Area = trMeans, StdErr = stdErrSingleFRRC, DF = dfSingleFRRC, CI = CISingleFRRC, row.names = NULL)
    colnames(ciAvgRdrEachTrtFRRC) <- c("Treatment", "Area", "StdErr", "DF", "CI Lower", "CI Upper")
    
    ssTFRRC <- array(0, dim = c(J))
    ssCFRRC <- array(0, dim = c(J))
    ssTCFRRC <- array(0, dim = c(J))
    for (j in 1:J) {
      for (i in 1:I) {
        ssTFRRC[j] <- ssTFRRC[j] + (mean(pseudoValues[i, j, ]) - mean(pseudoValues[, j, ]))^2
      }
      ssTFRRC[j] <- ssTFRRC[j] * K
      
      for (k in 1:K) {
        ssCFRRC[j] <- ssCFRRC[j] + (mean(pseudoValues[, j, k]) - mean(pseudoValues[, j, ]))^2
      }
      ssCFRRC[j] <- ssCFRRC[j] * I
      
      for (i in 1:I) {
        for (k in 1:K) {
          ssTCFRRC[j] <- ssTCFRRC[j] + (mean(pseudoValues[i, j, k]) - mean(pseudoValues[i, j, ]) - mean(pseudoValues[, j, k]) + mean(pseudoValues[, j, ]))^2
        }
      }
    }
    sourceArrayFRRC <- c("T", "C", "TC")
    dfArrayFRRC <- c(I - 1, K - 1, (I - 1) * (K - 1))
    ssArrayFRRC <- t(cbind(ssTFRRC, ssCFRRC, ssTCFRRC))
    ssTableFRRC <- data.frame(sourceArrayFRRC, dfArrayFRRC, ssArrayFRRC, row.names = NULL)
    colnames(ssTableFRRC) <- c("Source", "DF", readerID)
    
    msArrayFRRC <- ssArrayFRRC
    for (n in 1:3) msArrayFRRC[n, ] <- ssArrayFRRC[n, ]/dfArrayFRRC[n]
    msTableFRRC <- data.frame(sourceArrayFRRC, dfArrayFRRC, msArrayFRRC, row.names = NULL)
    colnames(msTableFRRC) <- c("Source", "DF", readerID)
    
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
    stdErrFRRC <- sqrt(2 * msArrayFRRC[3, ]/K)
    stdErrFRRC <- rep(stdErrFRRC, choose(I, 2))
    dim(stdErrFRRC) <- c(J, choose(I, 2))
    stdErrFRRC <- as.vector(t(stdErrFRRC))
    readerNames <- rep(readerID, choose(I, 2))
    dim(readerNames) <- c(J, choose(I, 2))
    readerNames <- as.vector(t(readerNames))
    trNames <- rep(diffTRName, J)
    dfReaderFRRC <- rep(K - 1, length(stdErrFRRC))
    CIReaderFRRC <- array(dim = c(length(stdErrFRRC), 2))
    tStat <- vector()
    tPr <- vector()
    for (n in 1:length(stdErrFRRC)) {
      tStat[n] <- diffTRMeansFRRC[n]/stdErrFRRC[n]
      tPr[n] <- 2 * pt(abs(tStat[n]), dfReaderFRRC[n], lower.tail = FALSE)
      CIReaderFRRC[n, ] <- sort(c(diffTRMeansFRRC[n] - qt(alpha/2, dfReaderFRRC[n]) * stdErrFRRC[n], diffTRMeansFRRC[n] + qt(alpha/2, dfReaderFRRC[n]) * stdErrFRRC[n]))
    }
    ciDiffTrtEachRdr <- data.frame(Reader = readerNames, Treatment = trNames, Estimate = diffTRMeansFRRC, StdErr = stdErrFRRC, DF = dfReaderFRRC, t = tStat, p = tPr, CI = CIReaderFRRC)
    colnames(ciDiffTrtEachRdr) <- c("Reader", "Treatment", "Estimate", "StdErr", "DF", "t", "Pr > t", "CI Lower", "CI Upper")
    if (option == "FRRC")
      return(list(fomArray = fomArray, anovaY = anovaY, anovaYi = msSingleTable, varComp = varComp, 
                  fFRRC = fFRRC, ddfFRRC = ddfFRRC, pFRRC = pFRRC, ciDiffTrtFRRC = ciDiffTrtFRRC, ciAvgRdrEachTrtFRRC = ciAvgRdrEachTrtFRRC, 
                  ssAnovaEachRdr = ssTableFRRC, msAnovaEachRdr = msTableFRRC, ciDiffTrtEachRdr = ciDiffTrtEachRdr))
  }
  
  if (option %in% c("RRFC", "ALL")) {
    # ************ RRFC ****************
    if (J > 1) {
      msDenRRFC <- msTR
      fRRFC <- msNum/msDenRRFC
      ddfRRFC <- ((I - 1) * (J - 1))
      pRRFC <- 1 - pf(fRRFC, I - 1, ddfRRFC)
      stdErrRRFC <- sqrt(2 * msDenRRFC/J/K)
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
        stdErrSingleRRFC[i] <- sqrt(msDenSingleRRFC[i]/J/K)
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
    if (option == "RRFC")
      return(list(fomArray = fomArray, anovaY = anovaY, anovaYi = msSingleTable, varComp = varComp, 
                  fRRFC = fRRFC, ddfRRFC = ddfRRFC, pRRFC = pRRFC, ciDiffTrtRRFC = ciDiffTrtRRFC, ciAvgRdrEachTrtRRFC = ciAvgRdrEachTrtRRFC))
  }
  
  return(list(fomArray = fomArray, anovaY = anovaY, anovaYi = msSingleTable, varComp = varComp, 
              fRRRC = fRRRC, ddfRRRC = ddfRRRC, pRRRC = pRRRC, ciDiffTrtRRRC = ciDiffTrtRRRC, ciAvgRdrEachTrtRRRC = ciAvgRdrEachTrtRRRC, 
              fFRRC = fFRRC, ddfFRRC = ddfFRRC, pFRRC = pFRRC, ciDiffTrtFRRC = ciDiffTrtFRRC, ciAvgRdrEachTrtFRRC = ciAvgRdrEachTrtFRRC, 
              ssAnovaEachRdr = ssTableFRRC, msAnovaEachRdr = msTableFRRC, ciDiffTrtEachRdr = ciDiffTrtEachRdr, 
              fRRFC = fRRFC, ddfRRFC = ddfRRFC, pRRFC = pRRFC, ciDiffTrtRRFC = ciDiffTrtRRFC, ciAvgRdrEachTrtRRFC = ciAvgRdrEachTrtRRFC))
} 
