print.CGP <-
function(x,...){
  
  cat("Call:\n")
  print(x$call)
  cat("\n Lambda:\n")
  print(x$lambda)
  cat("\n Theta:\n")
  print(x$theta)
  cat("\n Alpha:\n")
  print(x$alpha)
  cat("\n Bandwidth:\n")
  print(x$bandwidth)
  
}
