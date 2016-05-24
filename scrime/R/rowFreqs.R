`rowFreqs` <-
function(x, levels=1:3, divide.by.n=FALSE, affy=FALSE, includeNA=FALSE,
		useNN=c("not", "only", "also"), check=TRUE){
	tab <- rowTables(x, levels=levels, affy=affy, includeNA=includeNA, useNN=useNN,
		check=check)
	n <- if(divide.by.n) ncol(x) else rowSums(tab)
	tab/n
}

