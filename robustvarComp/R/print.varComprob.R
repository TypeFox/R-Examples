print.varComprob <- function(x, ...) {
  cat("\n Method: ", x$control$method, "\n\n")
  cat("Fixed effects: \n")
  print.default(x$beta)
  cat("Random effect variances: \n")
  print.default(x$eta)
  cat("Residual variance: \n")
  print.default(x$eta0)  
  cat("Value of the objective function at convergence: \n")
  print.default(x$min)
}

print.varComprob.compositeS <- print.varComprob.compositeTau <- print.varComprob.compositeMM <- print.varComprob.S <- print.varComprob.Tau <- print.varComprob.MM <- print.varComprob
