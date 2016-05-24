coreChecks <-
function(betaTruth,
											 X,
											 N,
											 etaTerms,
											 expandX,
											 betaNames)
{
	##
	if(ncol(as.matrix(X)) == 1)
    return("* 'X' should have at least two columns")
	
	##
	if(sum(X[,1] != 1) > 0)
  	return("* 'X' requires the first column to represent the intercept")
  
	##
	if(expandX != "none")
	{
		##
		if(expandX == "all")
			colIndex <- 1:ncol(X)
		if(expandX != "all")
			colIndex <- c(1:ncol(X))[is.element(colnames(X), expandX)]
		##
		for(i in 2:ncol(X))
		{
			if(is.element(i, colIndex))
			{
				prob1 <- sum(ceiling(X[,i]) != floor(X[,i]))
				prob2 <- (min(X[,i]) != 0)
				prob3 <- (max(X[,i]) != (length(unique(X[,i])) - 1))
				prob4 <- (max(X[,i]) == 0)
				if((prob1 + prob2 + prob3 + prob4) > 0)
					return("* check that each variable is consistent with the {0,1,2,...} coding convention")
			}
		}
	}
	
  ##
  if(length(N) != nrow(X))
    return("* incompatible dimensions of 'X' and 'N'")
	
	##
	Xtemp <- X
	if(!is.null(etaTerms))
	{
		##
		cat("\nWARNING: Make sure that the elements of 'etaTerms' are in the same order as those of 'betaTruth'\n\n")
		##
		if(sum(!is.element(etaTerms, colnames(X))) > 0)
			return("* elements of 'etaTerms' are not in the design matrix 'X'")
		Xtemp <- X[, is.element(colnames(X), etaTerms)]
	}
	
  ##
  if(expandX == "all")  p <- sum(unlist(lapply(apply(Xtemp, 2, unique), FUN=length)) - 1) + 1
	if(expandX != "all") p <- ncol(Xtemp)
	if(length(betaTruth) != p)
		return("* invalid dimension of 'betaTruth'")
	  
	##
	if(!is.null(betaNames))
	{
		if(length(betaTruth) != length(betaNames))
			return("* 'betaTruth' and 'betaNames' are not of the same length")
		if(!is.character(betaNames))
			return("* elements of 'betaNames' are not character")
	}

	##
	return("")
}
