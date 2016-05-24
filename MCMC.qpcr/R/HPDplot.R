HPDplot <-
function(model,factors,factors2=NULL,ylimits=NULL,hpdtype="w",inverse=F,jitter=0,plot=T,grid=T,zero=T,...){
	if (inverse) inv=-1 else inv=1
	res1=res2=c(rep(0,length(model$Sol[,1])))
	names1=names2=c()
	allnames=names(posterior.mode(model$Sol))
	genes=allnames[-grep(":",allnames)]
	genes=genes[grep("gene",genes)]
	genes=sub("gene","",genes)
	if(genes[1]=="NORM") { genes=genes[-1] } 
	ngenes=length(genes)
	first=1
	for (f in factors){		
		pattern=paste('^\\w+:',f,'$',sep="")
		f1=grep(pattern,allnames)
		res1=res1+model$Sol[,f1]
		if (first==1) {
			names1=allnames[f1]
			first=0
		} else {
			names1=paste(names1," + ",allnames[f1],sep="")
		}
	}
	first=1
	for (f in factors2){		
		pattern=paste('^\\w+:',f,'$',sep="")
		f1=grep(pattern,allnames)
		res2=res2+model$Sol[,f1]
		if (first==1) {
			names2=allnames[f1]
			first=0
		} else {
			names2=paste(names2," + ",allnames[f1],sep="")
		}
	}
	if (is.null(factors2)==FALSE){
		names.final=paste("(",names1,") - (",names2,")",sep="")
	} else {
		names.final=names1
	}
	stat=(inv*(res1-res2))/log(2)
	means=apply(stat,2,mean)
	hpds=HPDinterval(stat)
	if (plot==F) {
		modes=posterior.mode(stat)
		sds=apply(stat,2,sd)
		zs=means/sds
		ps=2*(1-pnorm(abs(zs)))
		pmc=mcmc.pval(stat)
		table.final=data.frame(cbind("mode"=modes,"mean"=means,hpds,"pval.z"=ps,"pval.mcmc"=pmc))
		row.names(table.final)=names.final
		return (table.final)
		stop
	}
	if(is.null(ylimits)) ylimits=c(min(hpds[,1])-0.5,max(hpds[,2])+0.5)
	plot(c(1:ngenes)+jitter,means,xlim=c(0.5,ngenes+0.5),xaxt="n",ylim=ylimits,xlab="",ylab="log2(fold change)",mgp=c(2.3,1,0),...)
	if (hpdtype=="l") {
		lines(hpds[,1],lty=2,...)
		lines(hpds[,2],lty=2,...)
	}	else {
		arrows(c(1:ngenes)+jitter,means,c(1:ngenes)+jitter,hpds[,1],angle=90,code=2,length=0.03,...)
		arrows(c(1:ngenes)+jitter,means,c(1:ngenes)+jitter,hpds[,2],angle=90,code=2,length=0.03,...)
		if (grid==TRUE){
			verts=seq(0.5,ngenes+0.5,1)
			abline(v=verts,lty=3,col="grey60")
		} 		
	} 		
	if (zero==TRUE) abline(h=0,lty=3,col="grey60")
	axis(side=1,labels=genes,at=c(1:ngenes),las=2)
}
