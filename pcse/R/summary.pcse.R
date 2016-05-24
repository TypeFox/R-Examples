# 
# default summary method for class "pcse"
# 



summary.pcse <- function(object, ...){
  
  K <- length(object$b)
  results.matrix <- matrix(NA, K, 4)
  results.matrix[,1] <- object$b
  results.matrix[,2] <- object$pcse
  results.matrix[,3] <- object$tstats
  results.matrix[,4] <- object$pval
  row.names(results.matrix) <- names(object$b)
  colnames(results.matrix) <- c("Estimate", "PCSE", "t value", "Pr(>|t|)")
  other.vals <- matrix(NA, 1, 3)
  other.vals[1,1] <- object$nobs
  other.vals[1,2] <- object$nmiss
  other.vals[1,3] <- object$df
  row.names(other.vals) <- ""
  colnames(other.vals) <- c("# Valid Obs", "# Missing Obs", "Df")

  cat("\n", "Results:", "\n", "\n")
  print(results.matrix)
  cat("\n", "---------------------------------------------", "\n", "\n")
  cat("# Valid Obs = ", other.vals[1,1], "; # Missing Obs = ", other.vals[1,2],
      "; Degrees of Freedom = ", other.vals[1,3], ".", "\n", sep="")
}
  
