#' Calculate number of cases for specified number of readers J to achieve the desired power for ROC studies.
#' 
#' Calculate required number of cases to achieve the desired power for specified number of readers J and DBM or OR variability parameters.
#' 
#' @param J The number of readers.
#' @param varYTR The DBM pseudovalue treatment-by-reader variance component.
#' @param varYTC The DBM pseudovalue treatment-by-case variance component.
#' @param varYEps The DBM pseudovalue error variance.
#' @param cov1 The covariance of the FOM resampling estimates for same reader and different modalities.
#' @param cov2 The covariance of the FOM resampling estimates for different readers and same modalities.
#' @param cov3 The covariance of the FOM resampling estimates for different readers and different modalities.
#' @param varEps The variance of the FOM resampling estimates for same reader and same modalities.
#' @param msTR Treatment-by-reader mean square of the FOM.
#' @param KStar The number of cases in the pilot study, only required when using OR variability parameters
#' @param alpha The significance level of the study, default value is 0.05.
#' @param effectSize The effect size to be used in the study, default value is 0.05.
#' @param desiredPower The desired statistical power, default value is 0.8.
#' @param randomOption It can be \code{"ALL"}, \code{"READERS"} or \code{"CASES"}, which indicate predictions for (1) random readers and random cases, (2) random readers only and 
#' (3) random cases only.
#' 
#' @return A list of two elements:
#' @return \item{K}{The minimum number of cases to just achieve the desired statistical power.}
#' @return \item{power}{The predicted statistical power.}
#' 
#' @details 
#' To calculate the sample size, either the group of DBM variance components (\code{varYTR}, \code{varYTC}, and \code{varYEps}) or 
#' OR covariance matrix elements and mean squares and number of cases in pilot study should be specified. 
#' If both of them are given, DBM variance components are used and OR values are ignored. Specifically, either numeric 
#' values of \code{varYTR}, \code{varYTC}, \code{varYEps} can be supplied, or the function call must explicitly state 
#' \code{cov1 = value1}, \code{cov2 = value2}, \code{cov3 = value3}, \code{varEps = value4}, \code{msTR = value5}, 
#' \code{KStar = value6}, as is standard usage in R.
#'
#' 
#' @references 
#' Hillis, S. L., Obuchowski, N. A., & Berbaum, K. S. (2011). Power Estimation for Multireader ROC Methods: An Updated and Unified Approach. Acad Radiol, 18, 129-142.
#' 
#' Hillis, S. L., Obuchowski, N. A, Schartz, K. M., & Berbaum, K. S. (2005). A comparison of the Dorfman-Berbaum-Metz and Obuchowski-Rockette methods 
#' for receiver operating characteristic (ROC) data. Statistics in Medicine, 24(10), 1579-607.
#' 
#' @examples
#' ## Following is an example of sample size calculation with DBM variance components.
#' retDbm <- DBMHAnalysis(data = rocData, fom = "Wilcoxon")
#' effectSize <- retDbm$ciDiffTrtRRRC$Estimate
#' varCompDBM <- retDbm$varComp
#' varYTR <- varCompDBM$varComp[3]
#' varYTC <- varCompDBM$varComp[4]
#' varYEps <- varCompDBM$varComp[6]
#' SampleSizeGivenJ(J = 6, varYTR = varYTR, varYTC = varYTC, varYEps = varYEps, 
#'                  effectSize =effectSize)
#'                      
#' ## Following is an example of sample size calculation with OR variance components.
#' retOR <- ORHAnalysis(data = rocData, fom = "Wilcoxon", covEstMethod = "Jackknife")
#' effectSize <- retOR$ciDiffTrtRRRC$Estimate
#' varCompOR <- retOR$varComp
#' cov1 <- varCompOR$varCov[3]
#' cov2 <- varCompOR$varCov[4]
#' cov3 <- varCompOR$varCov[5]
#' varEps <- varCompOR$varCov[6]
#' msTR <- retOR$msTR
#' KStar <- 114
#' SampleSizeGivenJ(J = 6, cov1 = cov1, cov2 = cov2, cov3 = cov3, varEps= varEps, 
#'                  msTR = msTR, KStar = KStar, effectSize =effectSize)
#'
#' \dontrun{ 
#' ## Following is an example of sample size calculation with DBM variance components, 
#' ## and scanning the number of readers
#' retDbm <- DBMHAnalysis(data = rocData, fom = "Wilcoxon")                     
#' effectSize <- retDbm$ciDiffTrtRRRC$Estimate
#' varYTR <- retDbm$varComp$varComp[3]
#' varYTC <- retDbm$varComp$varComp[4]
#' varYEps <- retDbm$varComp$varComp[6]
#' effectSize <- retDbm$ciDiffTrtRRRC$Estimate
#' for (J in 6:10) {
#'  ret <- SampleSizeGivenJ(J, varYTR, varYTC, varYEps, effectSize =effectSize) 
#'  message("# of readers = ", J, " estimated # of cases = ", ret$K, ", predicted power = ",
#'     signif(ret$power,3), "\n")
#' }
#' 
#' ## Following is an example of sample size calculation with OR variance components, 
#' ## using bootstrap to estimate variance components
#' retOR <- ORHAnalysis(data = rocData, fom = "Wilcoxon", covEstMethod = "Bootstrap")
#' effectSize <- retOR$ciDiffTrtRRRC$Estimate
#' varCompOR <- retOR$varComp
#' cov1 <- varCompOR$varCov[3]
#' cov2 <- varCompOR$varCov[4]
#' cov3 <- varCompOR$varCov[5]
#' varEps <- varCompOR$varCov[6]
#' msTR <- retOR$msTR
#' KStar <- length(rocData$NL[1,1,,1])
#' SampleSizeGivenJ(J = 6, cov1 = cov1, cov2 = cov2, cov3 = cov3, varEps= varEps, 
#'                  msTR = msTR, KStar = KStar, effectSize =effectSize)
#' }
#' 
#' @export
#' 
SampleSizeGivenJ <- function(J, varYTR, varYTC, varYEps, cov1, cov2, cov3, varEps, msTR, KStar, 
                             alpha = 0.05, effectSize = 0.05, desiredPower = 0.8, randomOption = "ALL") {
  if (((missing(varYTR)) || (missing(varYTC)) || (missing(varYEps)) || is.null(varYTR) || is.null(varYTC) || is.null(varYEps)
       || is.na(varYTR) || is.na(varYTC) || is.na(varYEps))) {
    if (((!missing(cov1)) && (!missing(cov2)) && (!missing(cov3)) && (!missing(varEps)) && (!missing(msTR)) && (!missing(KStar)))) {
      varYTR <- msTR - varEps + cov1 + max(cov2 - cov3, 0)
      varYTC <- max(cov2 - cov3, 0) * KStar
      varYRC <- max(cov1 - cov3, 0) * KStar
      varYEps <- (varEps - cov2) * KStar - varYRC
    } else {
      stop("Either DBM or OR variance components should be specified.")
    }
  }
  
  K <- 1
  power <- 0
  if (randomOption == "CASES") {
    if (varYTR != 0) 
      #warning("Readers effect is fixed. varYTR is reseted to 0.")
    varYTR <- 0
    while (power < desiredPower) {
      if (K > 2000) {
        break
      }
      K <- K + 1
      power <- PowerGivenJK(J, K, varYTR, varYTC, varYEps, alpha = alpha, effectSize = effectSize, randomOption = randomOption)
    }
  } else if (randomOption == "READERS") {
    if (varYTC != 0) 
      #warning("Cases effect is fixed. varYTC is reseted to 0.") 
    varYTC <- 0
    while (power < desiredPower) {
      if (K > 2000) {
        break
      }
      K <- K + 1
      power <- PowerGivenJK(J, K, varYTR, varYTC, varYEps, alpha = alpha, effectSize = effectSize, randomOption = randomOption)
    }
  } else {
    while (power <= desiredPower) {
      if (K > 2000) {
        break
      }
      K <- K + 1
      power <- PowerGivenJK(J, K, varYTR, varYTC, varYEps, alpha = alpha, effectSize = effectSize, randomOption = randomOption)
    }
  }
  return(list(K = K, power = power))
} 
