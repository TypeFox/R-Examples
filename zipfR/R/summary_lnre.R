summary.lnre <- function (object, ...)
{
  if (! inherits(object, "lnre")) stop("argument must belong to a subclass of 'spc'")

  object$util$print(object)             # name of model + parameters
  cat("Population size: S =", object$S, "\n")

  cat("Sampling method: ");             # type of sampling / approximate calculations
  if (object$multinomial) cat("multinomial") else cat("Poisson")
  cat(", ");
  if (object$exact) cat("with exact calculations.") else cat("approximations are allowed.")
  cat("\n")
  
  if ("bootstrap" %in% names(object)) cat("Bootstrapping data available for", length(object$bootstrap), "replicates\n")

  cat("\n")

  spc <- object$spc
  if (!is.null(spc)) { # estimated from observed spectrum -> show comparison & gof
    N <- N(spc)
    cat("Parameters estimated from sample of size N = ", N, ":\n", sep="")
    report <- cbind(c(V(spc), EV(object, N)),
                    rbind(Vm(spc, 1:5), EVm(object, 1:5, N)))
    colnames(report) <- c("V", "V1", "V2", "V3", "V4", "V5")
    rownames(report) <- c("   Observed:", "   Expected:")
    report <- as.data.frame(round(report, digits=2))
    report[[" "]] <- c("...", "...")
    print(report)
    cat("\n")
    
    cat("Goodness-of-fit (multivariate chi-squared test):\n")
    gof <- object$gof
    rownames(gof) <- "  "  # should produce nice alignment with other displays
    print(gof)
  } else {
    cat("(Model parameters have not been estimated from observed spectrum.)\n")
  }

  invisible(NULL)
}
