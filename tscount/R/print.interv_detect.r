print.interv_detect <- function(x, ...){
  if(!is.null(x$error_message)){
    warning(x$error_message)
    cat("\n")
  }else{
    cat("\n\tDetecting an intervention of given type at unknown time\n")
    cat("\nMaximum test statistic: ", x$test_statistic, " at time ", x$tau_max, sep="")
    if(!is.null(x$p_value)) cat(", p-value: ", x$p_value, sep="")
    cat("\n\n")
    if(!is.null(x$fit_interv)){
      cat("Fitted model with specified intervention at time ", x$tau_max, ":\n", sep="")
      print(x$fit_interv, ...)
    }
    invisible(x)
  }
}
