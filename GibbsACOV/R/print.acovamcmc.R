#' @param x A \code{acovamcmc} object
#' @param ... Ignored
#' @method print acovamcmc

print.acovamcmc <- function(x, ...){
  cat("Gibbs Sampler Iterations:",x$Iterations, "\n\n")
  cat("Convergence Status:",x$Convergence_Diag,"\n\n")
  cat("Gelman-Rubin Threshold for Convergence: <=",x$Gelman_Rubin_Threshold,"\n\n")
  cat(100*x$Credible_Interval_Coverage,"% Credible Intervals (lower, estimate, upper):","\n")
  print.default(format(x$Credible_Interval),print.gap=2L,quote=FALSE)
  cat("\n")
  cat("Run Time (total elapsed seconds):",x$Run_Time[3],"\n")
}
