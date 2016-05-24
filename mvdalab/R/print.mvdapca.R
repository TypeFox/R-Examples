print.mvdapca <- function (x, ...) {
  object <- x
  nobj <- nrow(x$scores)
  nvars <- ncol(x$Xdata)
  cat("Principal Component Analysis\n")
  cat("\nFit Summary: \nNumber of objects =", nobj, "\nNumber of Variables =", nvars)
  cat("\nPercent Variation Explained:\n")
  print(round(x$Percents.Explained, 3))
  cat("\nCross-Validation Results:\n")
  print(data.frame(MSEP = round(x$GVC, 3)))
  cat("\nEigenvalues:\n")
  print(round(diag(x$D), 3))
}
