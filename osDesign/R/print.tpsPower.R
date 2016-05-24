print.tpsPower <-
function(x, ...)
{
	##
  cat("Number of simulations, B:",x$B,"\n")
  ##
  cat("Phase I stratification variable(s):", x$strataNames, "\n")
	##
  if(is.null(x$NI))
    cat("Sample size at Phase I:", sum(x$N), "\n")
  else
  {
    cat("Sample size at Phase I for controls:", x$NI[1], "\n")
    cat("Sample size at Phase I for casess:", x$NI[2], "\n")
  }
  ##
  cat("\n'True' regession coefficients, betaTruth:")
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
