

#########################################################
within.NestedImputationList <- function (data, expr, ...){
   res <- data
   imp <- res$imputations
   CALL <- res$call
   NB <- length(imp)   
   NW <- length(imp[[1]])         
   for (ii in 1:NB){
      for (ww in 1:NW){
	   # this function is simply a copy of within.data.frame
		parent <- parent.frame()
		data <- imp[[ii]][[ww]]
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
		imp[[ii]][[ww]] <- data
					}
			  }									
	res$imputations <- imp		
	res <- NestedImputationList(imp)
	res$call <- CALL
	return(res)
}
###############################################################