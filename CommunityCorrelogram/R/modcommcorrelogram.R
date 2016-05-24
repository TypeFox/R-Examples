sig.gauss<-function(h,Ch,Cc,Cw){
	pred<-NA
	if(h<Cc) pred<-0
	else pred<-Ch*(1-exp(-Cw*(h-Cc)^2))
	return(pred)
}

sig.gauss.vec<-function(hvec,Ch,Cc,Cw){
	sapply(hvec,FUN=sig.gauss,Ch,Cc,Cw)
}

min.fun<-function(params,commcorr,alpha,pw){
	sum((commcorr[,4]-
		sig.gauss.vec(commcorr[,1],
			params[1],params[2],params[3]))^2*
		commcorr[,2]*(exp(-pw*(commcorr[,4]-alpha)^2)))
}


mod.commcorrelogram<-function(object,Ch=1,Cc=5,Cw=0.01,plot=T,alpha=0.05
	,alternative='one.tailed',pw=5,lgpos='topleft',...){
	commcorr<-na.omit(object@community.correlogram)
	if(alternative=='two.tailed'){
		commcorr[commcorr[,3]<0,4]<-1-commcorr[commcorr[,3]<0,4]
	}
	mlag<-max(commcorr$lag.distance)
	if(is.null(Cc)) Cc<-mlag/2
	temp.mod<-nlminb(c(1,2,0.03),min.fun,commcorr=commcorr
		,alpha=alpha,pw=pw,lower = 0, upper = Inf)
	pred<-data.frame(lag.distance=seq(0,mlag,0.05)
		,predicted.sig=sig.gauss.vec(seq(0,mlag,0.05)
		,temp.mod$par[1],temp.mod$par[2],temp.mod$par[3]))
	f<-function(h){
		alpha-sig.gauss.vec(h,temp.mod$par[1],temp.mod$par[2]
			,temp.mod$par[3])
	}
	f2<-function(h){
		(1-alpha)-sig.gauss.vec(h,temp.mod$par[1],temp.mod$par[2]
			,temp.mod$par[3])
	}
	range<-uniroot(f,c(0,mlag))$root
	rangeu<-NA
	if(alternative=='two.tailed'){
		rangeu<-uniroot(f2,c(0,mlag))$root
	}
	if(plot==TRUE){
		par(mfcol=c(2,1),mar=c(4,4,2,2))
		plot(commcorr[,c(1,3)]
			,pch=as.numeric(!(commcorr[,4]<alpha | commcorr[,4]>(1-alpha)))*2+19
			,xlab='lag distance',main='Significance Modeling',type='o',...)
		abline(h=0,lty=1,lwd=0.5)
		par(mar=c(5,4,1,2))
		plot(commcorr$lag.distance,commcorr[,4],xlab='lag distance',
			ylab='significance',pch=as.numeric(!(commcorr[,4]<alpha | commcorr[,4]>(1-alpha)))*2+19,ylim=c(0,1))
		lines(pred,col='blue')
		abline(h=alpha,lty=2,col='red')
		abline(v=range,col='red',lty=3)
		legend(lgpos,legend=c('empirical','model','alpha'
			,paste('range',round(range,digits=3),sep='=')),
			pch=c(19,-1,-1,-1),lty=c(-1,1,2,3),col=c('black','blue','red','red'),bg='white')}
	return(list(model.coefficients=temp.mod$par,empirical=commcorr,predicted=pred
		,range=range,outer.range=rangeu))
}
