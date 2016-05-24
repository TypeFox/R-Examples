print.interv_test <- function(x, ...){
  if(!is.null(x$error_message)){
    warning(x$error_message)
    cat("\n")
  }else{
    cat("\n\tScore test on intervention(s) of given type at given time\n")
    cat("\nChisq-Statistic: ", x$test_statistic, " on ", x$df, " degree(s) of freedom, p-value: ", x$p_value, "\n\n", sep="")
    if(!is.null(x$fit_interv)){
      cat("Fitted model with the specified intervention:\n")
      print(x$fit_interv, ...)
    }
  invisible(x)
  }
}
