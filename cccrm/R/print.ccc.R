print.ccc<-
function(x,...)
{
cat("CCC estimated by variance components: \n")
print(x$ccc[1:4])
cat("\n")
}
