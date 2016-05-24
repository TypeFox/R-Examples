
ssize.twoSamp<-function(delta,sigma,fdr=0.05,power=0.8,pi0=0.95,maxN=35,side="two-sided",cex.title=1.15,cex.legend=1){
	ds<-delta/sigma
	N<-maxN
	a<-fdr

	if(side=="two-sided"){
		TSfun2<-function(c){
			r<-a*(1-p)/((1-a)*p)
			dif<-abs((2*pt(q=-c,df=2*n-2)/(1-pt(q=c,df=2*n-2,ncp=t))-r))
			return(dif)
		}
	}

	if(side=="upper"){
		TSfun2<-function(c){
			r<-a*(1-p)/((1-a)*p)
			dif<-abs((pt(q=-c,df=2*n-2)/(1-pt(q=c,df=2*n-2,ncp=t))-r))
			return(dif)
		}
	}

	if(side=="lower"){
		TSfun2<-function(c){
			r<-a*(1-p)/((1-a)*p)
			dif<-abs((pt(q=-c,df=2*n-2)/pt(q=-c,df=2*n-2,ncp=t)-r))
			return(dif)
		}
	}
	
	pwr2<-NULL
	crit<-NULL
	ssize<-matrix(0,nrow=length(pi0),ncol=3)
	colnames(ssize)<-c("pi0", "ssize","power")
	up.start<-50

	for(i in 1:length(pi0)){
		p<-pi0[i]
		up<-up.start
		for(n in 2:N){
			ni<-n
			t<-ds/sqrt(2/ni)
			ci<-optimize(f=TSfun2,interval=c(0,up))$min
			up<-ci

			if(abs(ci-up.start)>=1){
				if(side=="two-sided"|side=="upper"){pwr.new<-1-pt(q=ci,df=2*ni-2,ncp=t)}
				if(side=="lower"){pwr.new<-pt(q=-ci,df=2*ni-2,ncp=t)}
			}
			if(abs(ci-up.start)<1){pwr.new<-0; ci<-NA}

			crit<-c(crit,ci)
			pwr2<-c(pwr2,pwr.new)

			if(pwr2[(i-1)*(N-1)+n-1]>=power & ssize[i,1]==0){	##finding first sample size with 
				ssize[i,]<-c(p,n,pwr2[(i-1)*(N-1)+n-1])		##power greater than desired power
			}
		}
		
	}

	ssize[,1]<-pi0
	if(sum(ssize==0)>0){warning("Desired power not achieved for at least one pi0")}
	ssize[ssize==0]<-NA

	pwrMatrix<-matrix(c(2:N,pwr2),ncol=length(pi0)+1,byrow=FALSE)
		
	for(i in 1:length(pi0)){
		if(i==1){
			plot(2:N,pwrMatrix[,i+1],col=i,xlim=c(0,N),ylim=c(0,1),xlab="",ylab="",pch=16)
			lines(2:N,pwrMatrix[,i+1],col=i,lty=i)
		}

		if(i!=1){
			points(2:N,pwrMatrix[,i+1],col=i,pch=16)
			lines(2:N,pwrMatrix[,i+1],col=i,lty=i)
		}

	}
	
	abline(h=power,lty=2,lwd=2)
	abline(v=0:N,h=0.1*(0:10),col="gray",lty=3)
	#title(xlab="Sample size (n)",ylab="Power",font.sub=4)
	title(xlab="Sample size (n)",ylab="Power")

	mtext(bquote(paste("Power vs. sample size with fdr=",.(round(fdr,4))," and ",
		Delta/sigma==.(round(ds,4)))),cex=cex.title,padj=-0.5)
	legend(x=N,y=0,xjust=1,yjust=0,col=1:i,pch=c(16,16,16),
		lty=1:length(pi0),legend=as.character(pi0),bg="white",
		title=expression(pi[0]),cex=cex.legend)

	pwrMatrix<-round(pwrMatrix,7)
	colnames(pwrMatrix)<-c("n",as.character(pi0))

	critMatrix<-matrix(c(2:N,crit),ncol=length(pi0)+1,byrow=FALSE)
	colnames(critMatrix)<-c("n",as.character(pi0))

	ret<-NULL
	ret$ssize<-ssize
	ret$power<-pwrMatrix
	ret$crit.vals<-critMatrix

	return(ret)
}





