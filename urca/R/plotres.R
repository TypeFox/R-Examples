##
## plotres
##
plotres <- function (x){
  if (!(class(x) == "ca.jo"))
    stop("\nObject is not of class 'ca.jo' \n")
  resids <- x@Z0 - x@Z1%*%t(x@GAMMA) - x@ZK%*%t(x@PI)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(1, 1))
  for (i in 1:x@P) {
    layout(matrix(c(1, 2, 1, 3), 2, 2))
    plot.ts(resids[, i], main = paste("Residuals of ", i, ". VAR regression", sep = ""), ylab = "", xlab = "")
    abline(h = 0, col = "red")
    acf(x@R0[, i], main = "Autocorrelations of Residuals")
    pacf(x@R0[, i], main = "Partial Autocorrelations of Residuals")
    if (interactive()){
      cat("\nType <Return> to continue: ")
      readline()
    }
  }
}
