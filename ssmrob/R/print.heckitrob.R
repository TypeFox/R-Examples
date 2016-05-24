print.heckitrob <-
function(x, ...)
{
  print(x$method)
  cat("selection stage coefficients: \n")
  print(x$stage1$coeff)
  cat("outcome coefficients: \n")
  print(x$stage2$coeff)
}
