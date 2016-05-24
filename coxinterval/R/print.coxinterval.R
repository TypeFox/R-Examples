print.coxinterval <- function(x, ...)
{
  cat("Call:\n")
  dput(x$call)
  digits <- max(3, getOption("digits") - 4)
  saved.digits <- options(digits = digits)
  on.exit(options(saved.digits))
  est <- x$coef
  se <- sqrt(diag(x$var))
  mat <- cbind(est, exp(est), se, est/se, 1 - pchisq((est/se)^2, 1))
  dimnames(mat) <- list(names(est), c("coef", "exp(coef)", "se(coef)", "z", "p"))
  cat("\n")
  if (x$p > 0) prmatrix(mat)
  else cat("Null model\n")
  cat("\n")
  if (!is.null(x$loglik)) {
    options(digits = ceiling(log10(x$n)) + digits)
    cat("Initial log-likelihood:", x$loglik[1], "\n")
    cat("Log-likelihood after", x$iter, "iterations:", x$loglik[x$iter + 1],
        "\n")
  }
  invisible(x)
}
