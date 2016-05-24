summary.ds <-
function(object, ...)
{
  x<-object
cat("\nInitial Data:\n ")
	print(x$IniDat)
cat("\n Command:\n")
 	print(x$Call)
cat("\n Analysis:\n")
	if(x$tipo=="O"){cat("Ordinary Dual Scaling\n")}
	if(x$tipo=="A"){cat("Force Classification DS\n")}
cat("\n Output Available on Demand:\n")
	summary(unclass(x))
}
