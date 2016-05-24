lensDetect <- function(x, crit=.05) {
  if(class(x) != "lensMod") {stop("x must be of class 'LensMod'")}
  val <- ifelse(x$'Validity p' < crit, 1, 0)
  uti <- ifelse(x$'Utilization p' < crit, 1, 0)
  dirMatch <- ifelse(x$'Cue Validities' > 0 & x$'Cue Utilizations' > 0, 1, ifelse(x$'Cue Validities' < 0 & x$'Cue Utilizations' < 0, 1, 0))
  flag <- ifelse(val==1 & uti==1 & dirMatch==0, 1, 0)
  tab <- list()
  for(i in 1:ncol(val)) {
    tab[[i]] <- table(factor(val[-1,i], levels=c(0,1))[flag[-1,i]==0], factor(uti[-1,i], levels=c(0,1))[flag[-1,i]==0])
  }
  sigFun <- function(L) {
    truePos <- L[2,2]
    trueNeg <- L[1,1]
    falsePos <- L[1,2]
    falseNeg <- L[2,1]
    Sens <- truePos / (truePos + falsePos)
    Spec <- trueNeg / (falsePos + trueNeg)
    Misint <- 1 - Spec
    Miss <- 1 - Sens
    likePos <- Sens / (1 - Spec)
    likeNeg <- (1 - Sens) / Spec
    posPP <- truePos / (truePos + falsePos)
    negPP <- trueNeg / (trueNeg + falseNeg)
    f <- 2 * ((posPP * negPP) / (posPP + negPP))
    Prev <- (truePos + falseNeg) / nrow(val)
    diagPow <- (falsePos + trueNeg) / nrow(val)
    corrClass <- (truePos + trueNeg) / nrow(val)
    falseClass <- 1 - corrClass
    oddsRatio <- (truePos * trueNeg) / (falsePos * falseNeg)
    Kappa <- cohen.kappa(L)$kappa
    sigOut <- c(truePos, trueNeg, falsePos, falseNeg, Sens, Spec, Misint, Miss, likePos, likeNeg, posPP, negPP, f, Prev, diagPow, corrClass, falseClass, oddsRatio, Kappa)
    sigOut
  }
  stats <- data.frame(matrix(unlist(lapply(tab, sigFun)), ncol=ncol(val), byrow=F), row.names=c("TruePos", "TrueNeg", "FalsePos", "FalseNeg", "Sensitivity", "Specificity", "Misinterpret", "MissRate", "PosLikelihood", "NegLikelihood", "PosPredPower", "NegPredPower", "F", "Prevalance", "DiagnosticPower", "CorrectClassRate", "FalseClassRate", "Odds-Ratio", "Kappa"))
  colnames(stats) <- colnames(x$'Cue Validities')
  stats
}