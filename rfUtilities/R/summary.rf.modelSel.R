#' @title Summarizing random forests model selection
#' @description Summarizing of the rf.modelSel function 
#'
#' @param object  Object of class rf.modelSel
#' @param ... Ignored
#'
#' @method summary rf.modelSel
#'
#' @export
summary.rf.modelSel <- function(object, ...) { 
  cat("Selected variables:", "\n")
    cat("\t", object$selvars, "\n")
    cat("\n")
  for(i in 1:length(object$parameters)) {
    cat("Variables in parameter set", i, "\n")
      cat("\t", object$parameters[[i]], "\n")
  	cat("\n")
  }
  cat("Variable importance for selected parameters:", "\n")
    cat("\n")
    print(object$importance)  
  	cat("\n")
  if( "rf.final" %in% ls(object) ) {
    cat("##################################", "\n")
    cat("Selected random forests model:", "\n")
    print(object$rf.final)
  }  
}
