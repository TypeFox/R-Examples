print.mlgarch <-
function(x, varma=FALSE, ...)
{
  pars <- coef.mlgarch(x, varma=varma)
  vcovmat <- vcov.mlgarch(x, varma=varma)
  out1 <- rbind(pars, sqrt(diag(vcovmat)))
  rownames(out1) <- c("Estimate:", "Std. Error:")
  out2 <- as.data.frame(matrix(NA,4,1))
  out2[1,1] <- as.character(round(logLik.mlgarch(x, varma=FALSE), digits=3))
  out2[2,1] <- as.character(round(logLik.mlgarch(x, varma=TRUE), digits=3))
  out2[3,1] <- as.character(round(x$aux$ynonzerorowsn, digits=0))
  out2[4,1] <- as.character(round(x$aux$yanyrowiszeron, digits=0))
  rownames(out2) <-   c("Log-likelihood (log-mgarch):",
    "Log-likelihood (varma):", "No. of obs. without zeros:",
    "No. of obs. with zeros:")
  colnames(out2) <- ""
  cat("\n")
  cat("Date:", x$date, "\n")
  cat("Method: Multivariate ML \n")
  cat("Message (nlminb):", x$message, "\n")
  cat("No. of observations:", x$aux$n, "\n")
  cat("Sample:", as.character(x$aux$y.index[1]),
    "to", as.character(x$aux$y.index[x$aux$n]), "\n")
  cat("\n")
  cat("Coefficients:\n")
  cat("\n")
  print(out1)
  print(out2)
  cat("\n")
}
