FDR=function(pvalues=NULL,alpha=0.95,dep=1){
 if (is.null(pvalues)) stop("No p-values entered")
	m=length(pvalues)
	if (dep<0) {const.m=sum(1/(1:m))} else {const.m=1}
	spvalues=sort(pvalues)
	FDR=(1-alpha)*(1:m)/(m*const.m)
	return(any(spvalues<FDR))
}
###
pvalue.FDR=function(pvalues=NULL,dep=1){
 if (is.null(pvalues)) stop("No p-values entered")
	m=length(pvalues)
	if (dep<0) {const.m=sum(1/(1:m))} else {const.m=1}
	spvalues=sort(pvalues)
	pv.FDR=min(spvalues*(m*const.m)/(1:m))
	return(pv.FDR)
}
###

