#' Pretty Printing of Hierarchical's Output
#'
#' Pretty printing of hierarchical's output.
#' 
#' @param x an object of class hierarchical
#' @param digits digits to round to
#' @param ... additional parameters
#' @usage \method{print}{hierarchical}(x, digits = 4, ...)
#' @return Invisible string of the printed object.
#' @export
#' @examples
#' 
#' # see ?hierarchical
#' 
#' 
print.hierarchical <- function(x, digits = 4, ...){
	
  ## argument checking and basic variable setting
  stopifnot(class(x) == "hierarchical")
  
  ## print call
  cat("Call:\n")  
  print(x$call)
  cat("\n")
  
  ## print method
  cat("Fitting method:\n")  
  if(x$method == "mcmc") cat("Metropolis-Hastings using Markov basis for proposals\n")
  if(x$method == "ipf") cat("Iterative proportional fitting (with stats::loglin)\n")
  cat("\n")
  
  ## print method
  cat("MCMC details:\n")  
  cat(paste0("N = ", x$iter, " samples (after thinning), burn in = ", x$burn, ", thinning = ", x$thin))
  cat("\n\n")  
  
  ## print statistics
  statMat <- cbind(
    "Stat" = round(x$statistic, digits), 
    "SE" = round(sapply(x$sampsStats, sd) / sqrt(x$iter), digits), 
    "P(>= stat)" = round(x$p.value, digits),
    "SE" = round(x$p.value.std.err, digits), 
    "mid-P-Value" = round(x$mid.p.value, digits)
  )
  statMat[1,1] <- statMat[1,2] <- ""
  statMat <- cbind(
    format(c("P(samp)", "Pearson X^2", "Likelihood G^2", 
      "Freeman-Tukey", "Cressie-Read", "Neyman X^2"), justify = "right"),
    statMat
  )
  statMat <- apply(statMat, 2, format, justify = "left")
  statMat <- rbind(
    c("Distance", "Stat", "SE", "p.value", "SE", "mid.p.value"),
    statMat
  )
  statMat <- apply(statMat, 2, format, justify = "right") 
  apply(statMat, 1, function(x){cat(x); cat("\n")})
  invisible()  
}



