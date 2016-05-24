print.summary.genfrail <- function(x, digits=4, ...) {
  sum.dat <- x
  
  if (!is.null(sum.dat$call$lambda_0)) {
    hazard.fn <- "lambda_0"
  } else if (!is.null(sum.dat$call$Lambda_0)) {
    hazard.fn <- "Lambda_0"
  } else if (!is.null(sum.dat$call$Lambda_0_inv)) {
    hazard.fn <- "Lambda_0_inv"
  } else {
    hazard.fn <- "(default) Lambda_0"
  }
  
  cat(
            "genfrail created     : ", strftime(sum.dat$created), "\n",
    sprintf("Observations         : %d", sum.dat$n.obs), "\n",
    sprintf("Clusters             : %d", sum.dat$n.clusters), "\n",
    sprintf("Avg. cluster size    : %.2f", sum.dat$mean.cluster), "\n",
    sprintf("Right censoring rate : %.2f", sum.dat$censor.rate), "\n",
    sprintf("Covariates           : %s(%s)", sum.dat$covar.distr, toString(sum.dat$covar.param)), "\n",
    sprintf("Coefficients         : %s", toString(format(round(sum.dat$beta, digits), nsmall=2))), "\n",
            "Frailty              : ",
    ifelse(sum.dat$frailty != "none", 
           sprintf("%s(%s)", sum.dat$frailty, toString(sum.dat$theta)), "none"), "\n",
            "Baseline hazard      : ", hazard.fn, "\n", 
            "                        = ", gsub("\n", "", toString(sum.dat$hazard)), "\n",
  sep="")
  
  invisible()
}