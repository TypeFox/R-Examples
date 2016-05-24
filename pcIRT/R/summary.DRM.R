#' @rdname drm
#'
#' @export
#'
#' @method summary DRM
#'
#' @param object object of class \code{DRM}
#' @param \dots \dots{}


summary.DRM <-
function(object, ...){
  
  cat("\n Call: ", deparse(object$call), "\n\n")
  
  cat("Function calls: \n")
  print(object$fun_calls)
  
  cat("Convergence: \t", ifelse(object$convergence==0, "convergence", "no convergence"), "\n\n")
  
  cat("Deviance: \t", deparse(round(object$logLikelihood*(-2),3)), "\n")
  cat("Number of Parameters: \t", length(object$estpar), "\n\n")
  
  cat("-----------------------------------------------------\n")
  
    parall <- rbind(cbind("item estimates"=object$itempar, "SE"=object$itempar_se, "SE low"= object$itempar_se_low, "SE up"=object$itempar_se_up),"distrpar"=c(object$distrpar, object$distrpar_se, object$distrpar_se_low, object$distrpar_se_up) )
  
  cat("Parameter estimates: \n")
  print(parall)  
  
}
