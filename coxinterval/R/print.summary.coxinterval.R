print.summary.coxinterval <-
  function(x, digits = max(3, getOption("digits") - 4),
           signif.stars = getOption("show.signif.stars"), ...)
{
  cat("Call:\n")
  dput(x$call)
  saved.digits <- options(digits = digits)
  on.exit(options(saved.digits))
  f <- function(fit, formula = TRUE) {
    if (formula) {
      cat("\n")
      cat("Formula:\n")
      dput(fit$formula)
    }
    if (is.null(fit$coef) & x$p > 0) cat("  Estimation failed\n")
    else {
      cat("\n")
      if (fit$p > 0) {
        printCoefmat(fit$coef, digits = digits, signif.stars = signif.stars,
                     dig.tst = max(1, digits - 1), ...)
        if (!is.null(fit$conf)) {
          cat("\n")
          prmatrix(fit$conf)
        }
      }
      else cat("Null model\n")
      cat("\n")
      cat("Based on ")
      cat("n =", fit$n)
      if (!is.null(fit$m))
        cat(" subjects contributing", fit$m, "observation times")
      if (length(fit$na.action)) {
        if (!is.null(fit$m)) cat("\n(")
        else cat(" (")
        cat(length(fit$na.action), "deleted due to missingness)\n")
      }
      else cat("\n")
    }
  }
  f(x, formula = FALSE)
  options(digits = ceiling(log10(x$n)) + digits)
  cat("\n")
  cat("Initial log-likelihood:", x$loglik[1], "\n")
  cat("Log-likelihood after", x$iter, "iterations:", x$loglik[x$iter + 1], "\n")
  options(digits = digits)
  cat("\n")
  prmatrix(x$censor.rate)
  if (!is.null(x$rcfit)) {
    cat("\n")
    if (is.null(x$censor))
      cat("Estimation from imputed data via timereg's cox.aalen function\n")
    else if (x$censor == "right")
      cat("Estimation from right-censored data via survival's coxph function\n")
    else
      cat("Estimation from imputed data via survival's coxph function\n")
    f(x$rcfit)
  }
  invisible()
}
