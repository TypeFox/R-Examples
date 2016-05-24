print.tpsSim <-
function(x, ...)
{
	##
  cat("Number of simulations, B:",x$B,"\n")
  ##
  if(!is.list(x$strata))
  {
  	if(max(x$strata) == 0)
  	{
	  	cat("Phase I stratification variable(s):\n")
	  	cat("\tAll combinations of\n")
	  	for(i in 1:length(x$strataNames)) cat("\t", i, ":", x$strataNames[i], "\n")
  	}
	  if(max(x$strata) > 0)
		  cat("Phase I stratification variable(s):", x$strataNames, "\n")
  }
  if(is.list(x$strata))
  {
  	cat("Phase I stratification variable(s):\n")
  	temp <- sort(unique(unlist(x$strata)))
	  for(i in 1:length(temp)) cat("\t", temp[i], ":", x$strataNames[i], "\n")
  }
  ##
  if(is.null(x$NI))
    cat("Sample size at Phase I:", sum(x$N), "\n")
  else
  {
    cat("Sample size at Phase I for controls:", x$NI[1], "\n")
    cat("Sample size at Phase I for casess:", x$NI[2], "\n")
  }
  ##
  if(!is.list(x$strata))
  {
	  if(max(x$strata) == 0)
  	  cat("Sample size at Phase II, nII=c(nII0, nII1):", x$nII, "\n")
  	if(max(x$strata) > 0)
  	{
	  	cat("Sample size at Phase II:")
			print(matrix(rbind(x$nII0, x$nII1), nrow=2, ncol=length(x$nII0), dimnames=list(c("  controls, nII0: ", "     cases, nII1: "), rep("", length(x$nII0)))))
  	}
  }
  if(is.list(x$strata))
  	  cat("Sample size at Phase II, nII=c(nII0, nII1):", x$nII, "\n")
  ##
  cat("Sample size for the case-control design, nCC:", x$nCC, "\n")
  ##
  cat("\n'True' regession coefficients, betaTruth:")
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
#  cat(paste(round((1-x$alpha)*100, 1), "%", sep=""),"coverage probability\n")
#  print(round(x$results$betaCP, digits=x$digits))
  cat("Power\n")
  print(round(x$results$betaPower, digits=x$digits))
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
