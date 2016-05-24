
print.clogitboost <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("Relative influence score:\n")
  print(x$infscore)
}