#' @param object A \code{acovamcmc} object
#' @param ... Ignored
#' @method summary acovamcmc

summary.acovamcmc <- function(object, ...){
  cat("Convergence Status:",object$Convergence_Diag,"\n\n")
  cat("Gelman-Rubin Threshold for Convergence: <=",object$Gelman_Rubin_Threshold,"\n\n")
  cat(100*object$Credible_Interval_Coverage,"% Credible Intervals (lower, estimate, upper):","\n")
  print.default(format(object$Credible_Interval),print.gap=2L,quote=FALSE)
}
