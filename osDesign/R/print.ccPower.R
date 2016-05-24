print.ccPower <-
function(x, ...)
{
	##
  cat("Number of simulations, B:",x$B,"\n")
  ##
  cat("'True' regession coefficients, betaTruth:")
  temp <- matrix(x$betaTruth, ncol=1)
	rownames(temp) <- paste("  ", colnames(x$betaPower))
	colnames(temp) <- ""
	print(temp)
	##
	cat("\nPower\n")
	print(round(x$betaPower, digits=x$digits))
	##
	if(max(x$failed) > 0)
	{
	  cat("\nNumber of failed repititions")
		print(x$failed)
	}
  ##
  invisible()
}
