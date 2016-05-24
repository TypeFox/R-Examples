print.fitfrail <- function(x, digits=max(options()$digits - 4, 3), ...) {
  fit <- x
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
  cat("Call: ")
  print(fit$call)
  
  cat("\n")
  tmp <- data.frame(Covariate=names(fit$beta), Coefficient=fit$beta)
  
  if (!is.null(fit$se)) {
    tmp$SE <- fit$se.beta
  }
  
  print(format(tmp, width=8, justify="right", nsmall=digits), print.gap=5, row.names=FALSE, right=TRUE)
  
  cat("\n")
  cat("Frailty distribution   ", 
      fit$frailty, 
      "(", toString(format(fit$theta, nsmall=digits)), "), ",
      "VAR of frailty variates = ", format(fit$frailty.variance, nsmall=digits), 
      "\n",
      sep="")
  
  if (!is.null(fit$se)) {
    cat("Frailty parameter SE   ", 
        toString(format(fit$se.theta, nsmall=digits)),
        "\n",
        sep="")
  }
  
  cat("Log-likelihood         ", 
      format(fit$loglik, nsmall=digits), 
      "\n", 
      sep="")
  cat("Converged (method)     ",
      fit$iter,
      " iterations, ",
      format(fit$fit.time, nsmall=2),
      ifelse(fit$fitmethod == "score", " (solved score equations)", " (maximized log-likelihood)"),
      "\n",
      sep="")
  
  invisible()
}