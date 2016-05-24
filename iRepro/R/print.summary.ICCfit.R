print.summary.ICCfit <-
function(x, ...){
  cat("Intraclass correlation for interval-censored data: maximum likelihood estimates\n\n")
  print(x$estimates)
  cat("\n")
  
  cat("---\n")
  print(x$loglikelihood)  
}

