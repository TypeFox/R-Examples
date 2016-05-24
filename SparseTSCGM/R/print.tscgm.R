 
print.tscgm <- function(x, ...){
  print(x$theta)
  print(x$gamma)
  cat("Minimum IC type value is =", x$min.ic, "\n")
  cat("optimal tuning parameter: lambda1 =", x$lam1.opt, "\n")
  cat("optimal tuning parameter: lambda2 =", x$lam2.opt, "\n")
  cat("Sparsity level of precision matrix: ", x$s.theta, "\n")
  cat("Sparsity level of autoregressive coefficients matrix: ", x$s.gamma, "\n")
}