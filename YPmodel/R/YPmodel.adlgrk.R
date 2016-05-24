YPmodel.adlgrk <-
function(data=c(), Estimate=c(), ...)
{
	if(is.null(data)){
		stop(paste(fun.errorMessage('DataSet')))
	}
	Data <- YPmodel.inputData(data=data)

	if(is.null(Estimate)){
		Estimate <- YPmodel.estimate(data=data)
	}

	beta <- Estimate$beta
	r <- Estimate$r
	#-----------------------------------------------------------------#
	## Step 3. Longrank test
	#-----------------------------------------------------------------#
	## loading data
	#-----------------------------------------------------------------#
	data8 <- fun.adlgrk(beta,r,Data)
	t <- data8$t
	ro <- data8$ro
	q <- data8$q
	pval <- data8$pval

	Adlgrk <- list(pval=pval)

	class(Adlgrk) <- "YPmodel.adlgrk"
	Adlgrk$call <- match.call()

	return(Adlgrk)

}
