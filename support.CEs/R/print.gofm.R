print.gofm <- function(x, digits = getOption("digits"), ...)
{
# Name : print.gofm
# Title: print() for S3 class "bws.count"
# Arguments:
#  x            An object of S3 class "gofm".
#  digits       The number of significant digits.
#  ...          Arguments passed to format().

# display various measures

  cat("\nRho-squared =", format(x$RHO2, digits = digits, ...), "\n")
  cat("Adjusted rho-squared =", format(x$AdjRHO2, digits = digits, ...), "\n")
  cat("Akaike information criterion (AIC) =", format(x$AIC, digits = digits, ...), "\n")
  cat("Bayesian information criterion (BIC) =", format(x$BIC, digits = digits, ...), "\n")
  cat("Number of coefficients =", x$K, "\n")
  cat("Log likelihood at start =", format(x$LL0, digits = digits, ...), "\n")
  cat("Log likelihood at convergence =", format(x$LLb, digits = digits, ...), "\n\n")

  invisible(x)
}

