print.dbmssEnvelope <-
function(x, ...) {
  
  einfo <- attr(x, "einfo")

  type <- ifelse(einfo$global, "Global", "Local")
  cat(paste(type, "critical envelopes obtained from", einfo$nsim, "simulations of", deparse(attr(x, "ylab"))))
  cat(paste(" under the null hypothesis:", einfo$H0, "\n"))

  if (!is.null(attr(x, "simfuns"))) 
    cat(paste("(All", einfo$nsim, "simulated function values are stored in attr(,", dQuote("simfuns"), ") )\n"))

  cat(paste("Significance level of Monte Carlo test:", einfo$Alpha, "\n"))
  cat(paste("Data:", einfo$Yname, "\n"))

  # Do not call NextMethod but print.fv directly (print.envelope is skipped)
  print.fv(x, ...)
}
