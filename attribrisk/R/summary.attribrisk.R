#' Summarizes the attribrisk object.
#'
#' This is a method for the function summary for object of the class "attribrisk".  
#' 
#' @param object A attribrisk output object of class "attribrisk". 
#' @param ... further arguments passed to or from other methods. 
#'
#' @return The attribrisk object invisible flag set to prevent printing.
#' @export

#' @seealso \code{\link{attribrisk}}

summary.attribrisk <- function(object, ...)
  {
    digits <- getOption("digits");
    if (hasArg(digits)) {
       digits <- list(...)$digits
    }
      
    print(summary(object=object$fit))
    cat('\n####################################################\n\n')

    if(!is.null(object$attribrisk.var)) {

      cat('Attributable risk:',signif(object$attribrisk,digits=digits), '\n')

      if (class(object$attribrisk.var)=='list' & ! ('boot' %in% names(object$attribrisk.var)) ) {
        for (i in 1:length(object$attribrisk.var)) {
          boot.arg <- deparse(object$call$var.config[[i+ 1]])
          if ("boot" %in% names(object$attribrisk.var[[i]])){
            cat("    The standard error of AR using ", boot.arg, " = ", round(sqrt(var(object$attribrisk.var[[i]]$boot$t)), digits=digits), '\n')
          } else {
            cat("    The standard error of AR using ", boot.arg, " = ", round(sqrt(object$attribrisk.var[[i]]), digits=digits), '\n')
          }
        }
      } else {
         boot.arg <- deparse(object$call$var.config)
         if ("boot" %in% names(object$attribrisk.var)){
           cat("    The standard error of AR using ", boot.arg, " = ", round(sqrt(var(object$attribrisk.var$boot$t)), digits=digits), '\n')
         } else {
           cat("    The standard error of AR using ", boot.arg, " = ", round(sqrt(object$attribrisk.var), digits=digits), '\n')
         }
      }

    } else {
      cat('Attributable risk:',signif(object$attribrisk,digits=digits),'\n')
    }
    
    invisible(object)
  }

