expandCatX <-
function(X, expandX="all")
{
	##
	if(expandX == "all")
		colIndex <- 1:ncol(X)
	if(expandX != "all")
		colIndex <- c(1:ncol(X))[is.element(colnames(X), expandX)]
	
	## Assumes a check has been performed to make sure the columns that are to be expanded
	## adhere to the {0,1,2,...} coding convention
	##
	value <- matrix(X[,1], ncol=1, dimnames=list(1:nrow(X), "Int"))
	##
	if(ncol(X) > 1)
	{
		n.lev <- unlist(lapply(apply(X, 2, unique), FUN=length))
		for(i in 2:ncol(X))
		{
			if(!is.element(i, colIndex))
			{
				value <- cbind(value, X[,i])
				colnames(value)[ncol(value)] <- colnames(X)[i]
			}
			else
			{
				for(j in 1:(n.lev[i]-1))
				{
					value <- cbind(value, as.numeric(X[,i] == j))
					if(n.lev[i] == 2) colnames(value)[ncol(value)] <- colnames(X)[i]
					if(n.lev[i] > 2) colnames(value)[ncol(value)] <- paste(colnames(X)[i], ".", j, sep="")
				}
			}
		}
	}
	##
	return(value)
}
