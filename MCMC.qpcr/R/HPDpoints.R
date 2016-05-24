HPDpoints <-
function(model,factors,factors2=NULL,ylimits=NULL,hpdtype="w",inverse=F,jitter=0,...){
	allnames=names(posterior.mode(model$Sol))
	genes=allnames[-grep(":",allnames)]
	genes=genes[grep("gene",genes)]
	genes=sub("gene","",genes)
	ngenes=length(genes)
	if (inverse) inv=-1 else inv=1
	res1=res2=c(rep(0,length(model$Sol[,1])))
	for (f in factors){		
		pattern=paste('^\\w+:',f,'$',sep="")
		f1=grep(pattern,allnames)
		res1=res1+model$Sol[,f1]
	}
	for (f in factors2){		
		pattern=paste('^\\w+:',f,'$',sep="")
		f1=grep(pattern,allnames)
		res2=res2+model$Sol[,f1]
	}
	means=apply((inv*(res1-res2))/log(2),2,mean)
	hpds=HPDinterval(inv*(res1-res2))
	hpds[,1]=hpds[,1]/log(2)
	hpds[,2]=hpds[,2]/log(2)
	points(c(1:ngenes)+jitter,means,...)
	if (hpdtype=="l") {
		lines(hpds[,1],lty=2,...)
		lines(hpds[,2],lty=2,...)
	}	else {
		arrows(c(1:ngenes)+jitter,means,c(1:ngenes)+jitter,hpds[,1],angle=90,code=2,length=0.03,...)
		arrows(c(1:ngenes)+jitter,means,c(1:ngenes)+jitter,hpds[,2],angle=90,code=2,length=0.03,...)
	} 		
}
