print.interv_multiple <- function(x, ...){
  if(!is.null(x$error_message)){
    warning(x$error_message)
    cat("\n")
  }else{
    cat("\n\tDetecting multiple intervention of unknown types at unknown times\n")
    cat("\nDetected intervention(s):\n", sep = "")
    if(nrow(x$interventions)==0){
      cat("No interventions detected\n")
    }else{
      print(x$interventions, ...)
    }
    cat("\nFitted model with all detected interventions:\n")
    print(x$fit_interv, ...)
  invisible(x)
  }
}
