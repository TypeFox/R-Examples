print.ccSim <-
function(x, ...)
{
	##
  cat("Number of simulations, B:",x$B,"\n")
  ##
  cat("'True' regession coefficients, betaTruth:")
  temp <- matrix(x$betaTruth, ncol=1)
	rownames(temp) <- paste("  ", colnames(x$results$betaPower))
	colnames(temp) <- ""
	print(temp)
  cat("\n")
	##
  cat("Mean percent bias\n")
  print(round(x$results$betaMeanPB, digits=x$digits))
  cat("\n")
	##
  cat(paste(round((1-x$alpha)*100, 1), "%", sep=""),"coverage probability\n")
  print(round(x$results$betaCP, digits=x$digits))
  cat("\n")
	##
  cat("Relative uncertainty\n")
  print(round(x$results$betaRU, digits=x$digits))
  cat("\n") 
	##
	if(max(x$failed) > 0)
	{
	  cat("Number of failed repititions")
		print(x$failed)
	}
  ##
  invisible()
}
