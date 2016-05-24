generate.roc <-
function(vals,pval=TRUE){
	sens = seq(0,1,by=.1)
	if(pval){
		pows = calculate.power(vals,alpha=rev(sens),comp=`<`)
	}else{
		pows = calculate.power(vals,alpha=sens,comp=`>`)
	}
	pows = aperm(pows,c(2,3,1))

	return(pows)
}
