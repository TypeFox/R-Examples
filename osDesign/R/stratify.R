stratify <-
function(X, strata=NULL)
{
	##
	if(is.null(strata))
		value <- 1:nrow(X)
	##
	if(!is.null(strata) & length(strata) == 1)
		value <- X[,strata]
	##
	if(!is.null(strata) & length(strata) > 1)
	{
		if(is.element(1, strata))
			strata <- sort(strata)[-1]
		##
		temp <- apply(as.matrix(X[,strata]), 2, unique)
		if(is.matrix(temp))
			n.lev <- apply(apply(as.matrix(X[,strata]), 2, unique), 2, FUN=length)
		if(is.list(temp))
			n.lev <- unlist(lapply(apply(as.matrix(X[,strata]), 2, unique), FUN=length))
		##
		base10 <- cumsum(ceiling(log10(n.lev))) - 1
		base10 <- matrix(10^base10, nrow=nrow(X), ncol=length(strata), byrow=TRUE)
		value  <- apply(X[,strata]*base10, 1, sum)
	}
	##
	temp   <- value
  Slvls  <- unique(sort(temp))
  for(k in 1:length(Slvls)) value[temp == Slvls[k]] <- k
	##
	return(value)
}
