print.MY <-
function(x, ...)
 {
  if (!inherits(x, "MY"))
     stop("argument 'x' must be of class 'MY'")
  cat("\n") 
  if (x$ypow == 1 & x$xpow == 1)
    cat("Estimated adjusted ", x$ymeasure, " of the response variable:\n\n", sep = "")
    else
    cat("Estimated adjusted ", x$ymeasure, " of the response variable in the ", x$space, " space:\n\n", sep = "")
  print(x$M, ...)
  cat("\n")
 }
