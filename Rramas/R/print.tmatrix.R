print.tmatrix <-
function(x,...){  
  cat("transition matrix:\n")
  class(x) <-"matrix"
  print(x)
  cat("\n\n")
 }

