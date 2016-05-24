print.summary.Bvs <- function(x,...){
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Inclusion Probabilities:\n")
  print(x$summary)
  cat("---\n")
  cat("Code: HPM stands for Highest posterior Probability Model and\n MPM for Median Probability Model.\n ")
  if(x$method=="gibbs"){
    cat("Results are estimates based on the visited models.\n")
  }
}
