#' Calculate statistical power given numbers of readers J and cases K for ROC studies
#' 
#' Calculate the statistical power with the given number of readers, number of cases and DBM or OR variances components.
#' 
#' @param J The number of readers to be used in the calculation.
#' @param K The number of cases to be used in the calculation.
#' @param varYTR The DBM variance component of treatment(modality)-by-reader interaction term.
#' @param varYTC The DBM variance component of treatment(modality)-by-case interaction term.
#' @param varYEps The variance component of DBM error term.
#' @param cov1 The OR covariances of the figure of merit estimates of same reader and different modalities.
#' @param cov2 The OR covariances of the figure of merit estimates of same reader and different modalities.
#' @param cov3 The OR covariances of the figure of merit estimates of same reader and different modalities.
#' @param varEps The variance component of OR error term.
#' @param msTR Treatment(modality)-by-reader mean square of the figure of merit.
#' @param KStar See \link{SampleSizeGivenJ}.
#' @param alpha The significantce level.
#' @param effectSize The effect size to be used in the calculation.
#' @param randomOption The random option. It can be \code{"ALL"}, \code{"READERS"} or \code{"CASES"}, which indicate predictions for (1) random readers and random cases, (2) random readers only and 
#' (3) random cases only.
#' 
#' @return The statistical power with given components and condition.
#' 
#' @details 
#' To calculate the statistical power, either the group of DBM variance components (\code{varYTR}, \code{varYTC}, and \code{varYEps}) or 
#' OR variance components (\code{cov1}, \code{cov2}, \code{cov3}, \code{varEps}, \code{msTR} and \code{KStar}) should be specified. 
#' If both of them are given, DBM variance components are used and OR variance components are ignored.
#' 
#' @examples
#' ## Following is an example of sample size calculation with DBM variance componements.
#' retDbm <- DBMHAnalysis(data = rocData, fom = "Wilcoxon")
#' effectSize <- retDbm$ciDiffTrtRRRC$Estimate
#' varCompDBM <- retDbm$varComp
#' varYTR <- varCompDBM$varComp[3]
#' varYTC <- varCompDBM$varComp[4]
#' varYEps <- varCompDBM$varComp[6]
#' PowerGivenJK(J = 6, K = 251, varYTR = varYTR, varYTC = varYTC, 
#'              varYEps = varYEps, effectSize = effectSize)
#'                      
#' ## Following is an example of sample size calculation with OR variance componements.
#' retOR <- ORHAnalysis(data = rocData, fom = "Wilcoxon", covEstMethod = "Jackknife")
#' effectSize <- retDbm$ciDiffTrtRRRC$Estimate
#' varCompOR <- retOR$varComp
#' cov1 <- varCompOR$varCov[3]
#' cov2 <- varCompOR$varCov[4]
#' cov3 <- varCompOR$varCov[5]
#' varEps <- varCompOR$varCov[6]
#' KStar <- length(rocData$NL[1,1,,1])
#' msTR <- retOR$msTR
#' PowerGivenJK(J = 6, K = 251, cov1 = cov1, cov2 = cov2, cov3 = cov3, 
#'              varEps = varEps, msTR = msTR, KStar = KStar, effectSize = effectSize)
#' 
#' @references 
#' Hillis, S. L., Obuchowski, N. A., & Berbaum, K. S. (2011). Power Estimation for Multireader ROC Methods: An Updated and Unified Approach. Acad Radiol, 18, 129-142.
#' 
#' Hillis, S. L., Obuchowski, N. a, Schartz, K. M., & Berbaum, K. S. (2005). A comparison of the Dorfman-Berbaum-Metz and Obuchowski-Rockette methods for receiver operating characteristic (ROC) data. 
#' Statistics in Medicine, 24(10), 1579-607.
#' 
#' @export
#' @importFrom stats qf pf
#' 
PowerGivenJK <- function(J, K, varYTR, varYTC, varYEps, cov1, cov2, cov3, varEps, msTR, KStar, alpha = 0.05, effectSize = 0.05, randomOption = "ALL") {
  if (((missing(varYTR)) || (missing(varYTC)) || (missing(varYEps)))) {
    if (((!missing(cov1)) && (!missing(cov2)) && (!missing(cov3)) && (!missing(varEps)) && (!missing(msTR)) && (!missing(KStar)))) {
      varYTR <- msTR - varEps + cov1 + max(cov2 - cov3, 0)
      varYTC <- max(cov2 - cov3, 0) * KStar
      varYRC <- max(cov1 - cov3, 0) * KStar
      varYEps <- (varEps - cov2) * KStar - varYRC
    } else {
      stop("Either DBM or (OR variance components, msTR and KStar) must be specified.")
    }
  } 
  fDen <- (varYTR + 1 / K * (varYEps + J * varYTC))
  delta <- ((effectSize)^2 * J/2) / fDen
  if (randomOption == "ALL") {
    ddfH <- fDen^2/((varYTR + 1 / K * varYEps)^2/(J - 1))
  } else if (randomOption == "READERS") {
    ddfH <- J - 1
  } else if (randomOption == "CASES") {
    ddfH <- K - 1
  }
  fvalue <- qf(1 - alpha, 1, ddfH)
  power <- pf(fvalue, 1, ddfH, ncp = delta, FALSE)
  return(power)
} 
