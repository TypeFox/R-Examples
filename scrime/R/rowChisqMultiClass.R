`rowChisqMultiClass` <-
function(..., listTables=NULL, add.pval=TRUE, sameNull=FALSE){
	if(is.null(listTables))
		listTables <- list(...)
	if(!is.list(listTables))
		stop("listTables must be a list (of matrices).")
	if(length(listTables)<2)
		stop("At least two matrices need to be specified.")	
	cn <- checkTableList(listTables)
	if("NA"%in%cn){
		warning("The column named NA is removed from all matrices.")
		ids <- which(cn=="NA")
		listTables <- lapply(listTables, function(x) x[,-ids])
	}
	n.cat <- length(listTables)
	mat.n <- sapply(listTables, rowSums)
	n.obs <- rowSums(mat.n)
	tab <- listTables[[1]]
	nom <- tab * tab / mat.n[,1]
	for(i in 2:n.cat){
		tab <- tab + listTables[[i]]
		tmp <- listTables[[i]] * listTables[[i]] / mat.n[,i]
		nom <- nom + tmp
	}
	if(any(tab==0)){
		if(sameNull)
			stop("All variables must show the same number of variables.")
		df <- rowSums(tab>0) - 1
		tab[tab==0] <- 1
	}
	else
		df <- ncol(tab) - 1
	df <- df * (n.cat - 1)
	stats <- rowSums(nom/tab) - 1
	stats <- n.obs * stats
	if(!add.pval)
		return(stats)
	pval <- pchisq(stats, df, lower.tail=FALSE)
	pval[df==0] <- 1
	structure(list(stats=stats, df=df, rawp=pval))	
}

