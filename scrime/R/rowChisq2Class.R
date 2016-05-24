`rowChisq2Class` <-
function(cases, controls, add.pval=TRUE, sameNull=FALSE){
	if(!is.matrix(cases) | !is.matrix(controls))
		stop("cases and controls must be matrices.")
	if(any(dim(cases)!=dim(controls)))
		stop("cases and controls have not the same dimensions.")
	rn <- rownames(cases)
	if(any(rn!=rownames(controls)))
		stop("The row names differ between cases and controls.")
	cn <- colnames(cases)
	if(any(colnames(controls)!=cn))
		stop("The column names differ between cases and controls.")
	if("NA" %in% cn){
		ids <- which(cn=="NA")
		cases <- cases[,-ids]
		controls <- controls[,-ids]
		warning("The column named NA is removed from both cases and controls.")
	}
	if(any(cases<0))
		stop("All values in cases must be non-negative integers.")
	if(any(controls<0))
		stop("All values in controls must be non-negative integers.")
	n.cases <- rowSums(cases)
	n.controls <- rowSums(controls)
	n.obs <- n.cases + n.controls
	tab <- cases + controls
	if(any(tab==0)){
		if(sameNull)
			stop("All variables must show the same number of levels.")
		df <- rowSums(tab>0) - 1
		tab[tab==0] <- 1
	}
	else
		df <- ncol(tab) - 1
	cases <- cases * cases
	cases <- cases / n.cases
	controls <- controls * controls
	controls <- controls / n.controls
	cases <- cases + controls
	cases <- rowSums(cases/tab)-1
	cases <- n.obs*cases
	if(!add.pval)
		return(cases)
	pval <- pchisq(cases, df, lower.tail=FALSE)
	pval[df==0] <- 1
	structure(list(stats=cases, df=df, rawp=pval))
}

