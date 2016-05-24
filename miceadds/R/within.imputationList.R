

#########################################################
within.imputationList <- function (data, expr, ...){
   res <- data
   imp <- res$imputations
   M <- length(imp)   
   for (ii in 1:M){
	   # this function is simply a copy of within.data.frame
		parent <- parent.frame()
		data <- imp[[ii]]
		e <- evalq(environment(), data, parent)
		eval(substitute(expr), e)
		l <- as.list(e)
		l <- l[!sapply(l, is.null)]
		nD <- length(del <- setdiff(names(data), (nl <- names(l))))
		data[nl] <- l
		if (nD) 
			data[del] <- if (nD == 1) 
				NULL
			else vector("list", nD)
		imp[[ii]] <- data
			}
	res$imputations <- imp
	return(res)
}
###############################################################