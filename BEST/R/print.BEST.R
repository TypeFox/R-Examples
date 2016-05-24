print.BEST <- function(x, digits=4, ...) {
  # Somewhat less quick and dirty print method for BEST objects.

  # Sanity checks:
  if(!inherits(x, "data.frame"))
    stop("x is not a valid BEST object")
  if(ncol(x) == 3 && all(colnames(x) == c("mu","nu","sigma"))) {
    oneGrp <- TRUE
  } else if (ncol(x) == 5 && all(colnames(x) == c("mu1", "mu2","nu","sigma1","sigma2"))) {
    oneGrp <- FALSE
  } else {
    stop("x is not a valid BEST object")
  }

  Rhat <- attr(x, "Rhat")
  n.eff <- attr(x, "n.eff")
  doPriorsOnly <- attr(x, "doPriorsOnly")

  toPrint <- cbind(
    mean = colMeans(x),
    sd = apply(x, 2, sd),
    median = apply(x, 2, median), 
    t(hdi(x)))
  colnames(toPrint)[4:5] <- c("HDIlo", "HDIup")
  if(!is.null(Rhat))
    toPrint <- cbind(toPrint, Rhat = Rhat)
  if(!is.null(n.eff))
    toPrint <- cbind(toPrint, n.eff = round(n.eff))

  if(!is.null(doPriorsOnly) && doPriorsOnly) {
    cat("MCMC fit results for BEST: PRIORS ONLY!\n")
  } else {
    cat("MCMC fit results for BEST analysis:\n")
  }
  cat(nrow(x), "simulations saved.\n")
  print(toPrint, digits = digits)
  cat("\n'HDIlo' and 'HDIup' are the limits of a 95% HDI credible interval.\n")
  if(!is.null(Rhat))
    cat("'Rhat' is the potential scale reduction factor (at convergence, Rhat=1).\n")
  if(!is.null(n.eff))
    cat("'n.eff' is a crude measure of effective sample size.\n")

}
