print.loop2rsummary <- function (x,...) {
  cat("Summary Call:\n")
  print(x$summarycall)
  cat("Call for Original Fit:\n")
  print(x$call)
    cat("\nBootstrapped Estimates:\n")
    print(x$values[-c(1,2),c("Boot.Estimate","Bias","Std.Error","B.q0.025","B.q0.975")],digits=4)
 
}
