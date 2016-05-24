startingModel<-function(otree, odata, shifts=NULL, brownian = FALSE)
{
n<-otree@nterm;nt<-dim(odata)[2]
if(brownian){
	fit1<-apply(odata,2,function(x)brown(x,otree))
	LnL1<-sum(sapply(fit1,function(x)summary(x)$loglik))
	aic1<-getAIC(LnL1,2*nt,n*nt,TRUE)
	result<-list(list(fit=fit1,all_aic=NULL,aic=aic1[1], savedshifts = NULL, n_regimes=c(k=0,kprime=0,deltak=0,c=0,kprime_conv=0,kprime_nonconv=0)))
}else{
	if(!is.null(shifts)){
		if("1"%in%names(shifts)|"a"%in%shifts){
			if(names(shifts)[1]!="1"|shifts[1]!="a") stop("the first element of 'shifts' is automatically set to regime 'a' at node '1' - do not use these in 'shifts' except as the first element")
		if("1"%in%names(shifts))shifts<-shifts[-which(names(shifts)=="1")]
		if("a"%in%shifts)shifts<-shifts[-which(shifts=="a")]
		}
	}
	shifts1<-c("1"="a",shifts)
	regs1<-repaint(otree,shifts1)
	k<-length(shifts1);kk<-length(unique(shifts1))

	fit1<-list()
	rootage<-max(otree@times)
	for(i in 1:nt){
		initsigmasq<-var(odata[,i],na.rm=T)/rootage
		initalpha<-2/rootage
		fit1[[i]]<-hansen(odata[,i,drop=F],otree,regimes=regs1,sqrt.alpha=sqrt(initalpha),sigma=sqrt(initsigmasq))
	}
	names(fit1)<-names(odata)
	LnL1<-sum(sapply(fit1,function(x)summary(x)$loglik))
	aic1<-getAIC(LnL1,k+nt*(2+kk),n*nt,TRUE)
	xx<-summary(factor(shifts1))
	n_regimes<-c(k=k,kprime=kk,deltak=k-kk,c=sum(xx[xx>1]),kprime_conv=sum(xx>1),kprime_nonconv=sum(xx==1))
	result<-list(list(fit=fit1,all_aic=NULL,aic=aic1[1],savedshifts=shifts1,n_regimes=n_regimes))
	}
return(result)
}
