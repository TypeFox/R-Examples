maximizeOK <- function(fitobject=blackbox.getOption("fitobject"),
                       cleanResu="") {
  fittedNames <- blackbox.getOption("fittedNames")
  FONKgLow <- blackbox.getOption("FONKgLow")
  FONKgUp <- blackbox.getOption("FONKgUp")
  message.redef("\n*** Maximization ***")
  lower <- FONKgLow[fittedNames];upper <- FONKgUp[fittedNames]
  ## as.data.frame(xg[which.max(z), , drop=FALSE]) to keep names as 'names', (vs 'colnames' in matrix )
  initptinfK <- as.data.frame(fitobject$x[which.max(fitobject$fitted.values), , drop=FALSE])
  rosglobal <- findglobalMLE(initptinfK=initptinfK)
  blackbox.options(rosglobal=rosglobal) ## global
  if ("IBD" %in% blackbox.getOption("DemographicModel")) {
    range2Ns2 <- range(gridfn("latt2Ns2"))
    disttolow <- (rosglobal$latt2Ns2-min(range2Ns2))/(max(range2Ns2)-min(range2Ns2))
    if ( ! is.character( (extrascale <- blackbox.getOption("extraScale"))["latt2Ns2"] ) ## tested value may be NA or NULL dependeing on vector=NULL or not.
         && ! is.nan(disttolow) && disttolow<0.02) {
      .blackbox.data$options$extraScale <- c(extrascale, latt2Ns2="logscale")
    }
  }
  returncode <- rosglobal$convergence
  tmp <- rosglobal$edgelevel
  localst <- paste(blackbox.getOption("dataFile"), "(primary)", sep="") ## 26/04/11
  if (tmp>0) returncode <- returncode+tmp/(10^ceiling(log(tmp, 10))) ## second summand goes in decimal part of returcode
  writeoutput(localst, returncode=returncode, NA, NA, blackbox.getOption("CovFnParam")["smoothness"])
  message.redef("Primary likelihood estimates, and predicted Log(likelihood):")
  message.redef(prettyNamedUserValues(c(rosglobal$canonVP, "ln(L)"=blackbox.getOption("rosglobal")$value), extradigits=2)) ## canon is not in logscale...
  ## deriving an estimate of the number of points above the CI threshold (after maxim!)
  CIdlr <- -qchisq(1-blackbox.getOption("CIlevel"), 1)/2
  ptNbrforCI <- length(which(fitobject$fitted.values >
                               blackbox.getOption("rosglobal")$value+CIdlr))
  if (ptNbrforCI<(5^blackbox.getOption("fittedparamnbr"))) {
    if (ptNbrforCI==0) {
      warningue <- paste("(!) No computed point has a predicted likelihood above  the one-dimensional CI threshold.")
    } else warningue <- paste("(!) Only ", ptNbrforCI, " points have a predicted likelihood\n     in the upper ", signif(-CIdlr, 4), " [ln(L) units] range.", sep="")
    message.redef(warningue)
    write(warningue, file=cleanResu)
    warningue <- paste("    (this threshold corresponds to the ", signif(1-pchisq(-CIdlr*2, 1), 4), " chi-square threshold with 1 df);", sep="")
    message.redef(warningue)
    write(warningue, file=cleanResu)
    warningue <- "    It is advised to compute more points in order to obtain good CIs."
    message.redef(warningue)
    write(warningue, file=cleanResu)
  }
  blackbox.options(ptNbrforCI=ptNbrforCI)
  invisible(rosglobal)
}
