
#' @export 
#' 

print.MedianPolishM <-
function(x,...)
{
cat("Residuals:\n")
print(x$residuals)
cat("Overall:\n")
print(x$overall)
cat("\n")
cat("Effects:\n")
print(x$effects)
cat("Iter:\n")
print(x$iter)
}
