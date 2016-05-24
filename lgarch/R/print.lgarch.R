print.lgarch <-
function(x, arma=FALSE, ...)
{
  #out1:
  pars <- coef.lgarch(x, arma=arma)
  vcovmat <- vcov.lgarch(x, arma=arma)
  out1 <- rbind(pars, sqrt(diag(vcovmat)))
  rownames(out1) <- c("Estimate:", "Std. Error:")

  #out2:
  out2 <- as.data.frame(matrix(NA,4,1))
  out2[1,1] <- as.character(round(logLik.lgarch(x, arma=FALSE), digits=3))
  out2[2,1] <- as.character(round(logLik.lgarch(x, arma=TRUE), digits=3))
  out2[3,1] <- as.character(round(rss(x), digits=3))
  out2[4,1] <- as.character(round(x$aux$ynonzeron, digits=0))
  out2[5,1] <- as.character(round(x$aux$yzeron, digits=0))
  rownames(out2) <-   c("Log-likelihood (log-garch):",
    "Log-likelihood (arma):", "Sum Squared Resids. (arma):",
    "No. of non-zeros:", "No. of zeros:")
  colnames(out2) <- ""

  #header:
  cat("\n")
  cat("Date:", x$date, "\n")
  if(x$aux$method=="ls"){ methodmssg <- "LS" }
  if(x$aux$method=="ml"){ methodmssg <- "ML" }
  if(x$aux$method=="cex2"){ methodmssg <- "CEX2" }
  if(x$aux$mean.correction){
    methodmssg <- paste(methodmssg,
      " (w/mean-correction)", sep="")
  }
  cat("Method:", methodmssg, "\n")
  cat("Message (nlminb):", x$message, "\n")
  cat("No. of observations:", x$aux$n, "\n")
  cat("Sample:", as.character(x$aux$y.index[1]),
    "to", as.character(x$aux$y.index[x$aux$n]), "\n")
  cat("\n")
  cat("Coefficients:\n")
  print(out1)
  print(out2)
  cat("\n")
}
