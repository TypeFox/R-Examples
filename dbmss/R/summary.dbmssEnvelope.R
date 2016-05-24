summary.dbmssEnvelope <-
function(object, ...) {
  
  einfo <- attr(object, "einfo")

  type <- ifelse(einfo$global, "Global", "Local")
  cat(paste(type, "critical envelopes obtained from", einfo$nsim, "simulations of", deparse(attr(object, "ylab"))))
  cat(paste(" under the null hypothesis:", einfo$H0, "\n"))

  if (!is.null(attr(object, "simfuns"))) 
    cat(paste("(All", einfo$nsim, "simulated function values are stored in attr(,", dQuote("simfuns"), ") )\n"))

  cat(paste("Significance level of Monte Carlo test:", einfo$Alpha, "\n"))
  cat(paste("Data:", einfo$Yname, "\n"))

  return(invisible(NULL))
}