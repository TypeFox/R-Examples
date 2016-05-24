`print.SII` <-
function (x, digits=3, ...) 
{
  cat("\n")
  cat("SII:", round(x$sii, digits), "\n")
  cat("\n")  
}
