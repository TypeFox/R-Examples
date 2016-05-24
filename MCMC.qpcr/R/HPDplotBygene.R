
HPDplotBygene <-
function(model,gene,conditions,pval="mcmc",newplot=T,ylimits=NULL,inverse=F,jitter=0,plot=T,yscale="log2",interval="ci",grid=F,zero=F,...){
	if (inverse) inv=-1 else inv=1
	res1=res2=c(rep(0,length(model$Sol[,1])))
	names1=names2=c()
	allnames=names(posterior.mode(model$Sol))
	intercept=paste('gene',gene,sep="")	

	if (yscale=="log2") {
		model$Sol=model$Sol/log(2)
		Ylab="log2(abundance)"
	} else {
		if (yscale=="log10") {
			model$Sol=model$Sol/log(10)
			Ylab="log10(abundance)"
		}	else {
				Ylab="ln(abundance)" 
		}
	}
	
	cstats=c();means=c();hpds=c();pz=c()
	for (c in 1:length(names(conditions))) {
		co=names(conditions)[c]
		factors=conditions[[co]]$factors
		factors2=conditions[[co]]$factors2	
		res1=res2=c(rep(0,length(model$Sol[,1])))
		for (f in factors){		
			if (f==0) { 
				f1=intercept
			} else {	
				pattern=paste('gene',gene,':',f,"$",sep="")
				f1=grep(pattern,allnames)
			}
			res1=res1+model$Sol[,f1]
		}
		for (f in factors2){		
			if (f==0) { 
				f2=intercept
			} else {	
				pattern=paste('gene',gene,':',f,sep="")
				f2=grep(pattern, allnames)
			}
			res2=res2+model$Sol[,f1]
		}
		stat=(inv*(res1-res2))
		means=append(means,mean(stat))
		if (interval=="sd") hpds=rbind(hpds,c(mean(stat)-sd(stat),mean(stat)+sd(stat)),deparse.level=0)
		if (interval=="ci") hpds=rbind(hpds,HPDinterval(stat),deparse.level=0)
		row.names(hpds)=paste(row.names(hpds),co)
		ress=data.frame(cbind("mean"=means,"lo"=hpds[,1],"up"=hpds[,2]),deparse.level=0)
#print(ress)
		if (yscale=="proportion") {
			ress[ress<0]=0
			ress=sin(ress)^2
			Ylab="proportion" 
		}
		
#		sds=append(sds,sd(stat))
		cstats=cbind(cstats,stat,deparse.level=0)
		names(cstats)[c]=co
	}
	cstats=data.frame(cstats)
	# toplot=data.frame(cbind("mean"=means,"upper"=hpds[,1],"lower"=hpds[,2]))
	# toplot$condition=names(conditions)
	tukey=data.frame(diag(length(conditions))*0)
	names(tukey)=names(conditions)
	row.names(tukey)=names(conditions)
	meant=tukey
	for (cc in 1:(length(names(tukey))-1)) {
		for (cc2 in (cc+1):length(names(tukey))){
			stat=cstats[,cc2]-cstats[,cc]			
			zs=mean(stat)/sd(stat)
			if (pval=="z") pv=2*(1-pnorm(abs(zs)))	
			if (pval=="mcmc") pv=mcmc.pval(stat)	
			tukey[cc,cc2]=pv
			meant[cc,cc2]=mean(stat)
		}
	}
	if (plot==F) {
		return(list(mean.pairwise.differences=meant,pvalues=tukey,min=min(ress$lo),max=max(ress$up)))
		stop
	}
	marg=0.25
	if (yscale=="proportion") { marg=0.05 }
	if(is.null(ylimits)) ylimits=c(min(ress[,2])-marg,max(ress[,3])+marg)
	if (newplot==T) 	plot(c(1:length(conditions))+jitter,ress[,1],xlim=c(0.5,length(conditions)+0.5),xaxt="n",ylim=ylimits,xlab="",ylab=Ylab,mgp=c(2.3,1,0),...)
	if (newplot==F) points(c(1:length(conditions))+jitter,ress[,1],...)
	lines(c(1:length(conditions))+jitter,ress[,1],type="l",...)
	arrows(c(1:length(conditions))+jitter,ress[,1],c(1:length(conditions))+jitter,ress[,2],angle=90,code=2,length=0.03,...)
	arrows(c(1:length(conditions))+jitter,ress[,1],c(1:length(conditions))+jitter,ress[,3],angle=90,code=2,length=0.03,...)
		if (grid==TRUE){
			verts=seq(0.5,length(conditions)+0.5,1)
			abline(v=verts,lty=3,col="grey60")
		} 		 		
	if(newplot==T) axis(side=1,labels=names(conditions),at=c(1:length(conditions)),las=2)
	return(list(mean.pairwise.differences=meant,pvalues=tukey))
}
