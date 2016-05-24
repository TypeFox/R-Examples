
ssize.oneSampVary<-function(deltaMean,deltaSE,a,b,fdr=0.05,power=0.8,pi0=0.95,maxN=35,side="two-sided",cex.title=1.15,cex.legend=1){

	getAvgTcdf_varySigma<-function(c,a,b,deltaMean,deltaSE,n){
		sigmaFun<-function(rho){
			avgtcf<-(pt(q=c/sqrt(rho*deltaSE^2*n+1),df=n-1,
					ncp=deltaMean/sqrt((deltaSE^2+1/(rho*n))))
					*(1/(gamma(a)*(1/b)^a))*rho^(a-1)*exp((-1)*rho*b))
		}
		
		sigmaInt<-integrate(sigmaFun,0,Inf,abs.tol=1e-10)
		return(sigmaInt$value)
	}



	delMean<-deltaMean
	delSig<-deltaSE
	a<-a
	b<-b
	N<-maxN

	if(side=="two-sided"){
		OSVary<-function(c,fdr,p,n,a,b,dM,dS){
			r<-fdr*(1-p)/((1-fdr)*p)
			dif<-abs((2*pt(q=-c,df=n-1)/(1-getAvgTcdf_varySigma(c,a,b,dM,dS,n)
					+getAvgTcdf_varySigma(-c,a,b,dM,dS,n))-r))			
		}		
	}
	
	if(side=="upper"){
		OSVary<-function(c,fdr,p,n,a,b,dM,dS){
			r<-fdr*(1-p)/((1-fdr)*p)
			dif<-abs((pt(q=-c,df=n-1)/(1-getAvgTcdf_varySigma(c,a,b,dM,dS,n))-r))		
		}		
	}

	if(side=="lower"){
		OSVary<-function(c,fdr,p,n,a,b,dM,dS){
			r<-fdr*(1-p)/((1-fdr)*p)
			dif<-abs((pt(q=-c,df=n-1)/getAvgTcdf_varySigma(-c,a,b,dM,dS,n)-r))			
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
		for(n in 3:N){
			ci<-optimize(f=OSVary,interval=c(0,up),fdr=fdr,p=p,n=n,a=a,b=b,dM=delMean,dS=delSig)$min
			up<-ci

			if(abs(ci-up.start)>=1){
				if(side=="two-sided"){pwr.new<-(1-getAvgTcdf_varySigma(ci,a,b,delMean,delSig,n)
					+getAvgTcdf_varySigma(-ci,a,b,delMean,delSig,n))}
				if(side=="upper"){pwr.new<-1-getAvgTcdf_varySigma(ci,a,b,delMean,delSig,n)}
				if(side=="lower"){pwr.new<-getAvgTcdf_varySigma(-ci,a,b,delMean,delSig,n)}
			}
			if(abs(ci-up.start)<1){pwr.new<-0; ci<-NA}

			crit<-c(crit,ci)
			pwr2<-c(pwr2,pwr.new)

			if(pwr2[(i-1)*(N-2)+n-2]>=power & ssize[i,1]==0){	##finding first sample size with 
				ssize[i,]<-c(p,n,pwr2[(i-1)*(N-2)+n-2])		##power greater than desired power
			}
		}
	}

	ssize[,1]<-pi0
	if(sum(ssize==0)>0){warning("Desired power not achieved for at least one pi0")}
	ssize[ssize==0]<-NA

	pwrMatrix<-matrix(c(3:N,pwr2),ncol=length(pi0)+1,byrow=FALSE)
		
	for(i in 1:length(pi0)){
		if(i==1){
			plot(3:N,pwrMatrix[,i+1],col=i,xlim=c(0,N),ylim=c(0,1),xlab="",ylab="",pch=16)
			lines(3:N,pwrMatrix[,i+1],col=i,lty=i)
		}

		if(i!=1){
			points(3:N,pwrMatrix[,i+1],col=i,pch=16)
			lines(3:N,pwrMatrix[,i+1],col=i,lty=i)
		}

	}
	
	abline(h=power,lty=2,lwd=2)
	abline(v=0:N,h=0.1*(0:10),col="gray",lty=3)
	title(xlab="Sample size (n)",ylab="Power")
	mtext(bquote(paste("Average power vs. sample size with fdr=",.(fdr),",")),
		cex=cex.title,padj=-1.85)
	mtext(bquote(paste(Delta[g],"~N(",.(round(deltaMean,4)),",",.(round(deltaSE,4)),") and ",
		sigma[g]^2,"~IG(",.(round(a,4)),",",.(round(b,4)),")")),cex=cex.title,padj=-0.1)
	
	legend(x=N,y=0,xjust=1,yjust=0,col=1:i,pch=c(16,16,16),lty=1:length(pi0),
		legend=as.character(pi0),bg="white",title=expression(pi[0]),cex=cex.legend)

	pwrMatrix<-round(pwrMatrix,7)
	colnames(pwrMatrix)<-c("n",as.character(pi0))

	critMatrix<-matrix(c(3:N,crit),ncol=length(pi0)+1,byrow=FALSE)
	colnames(critMatrix)<-c("n",as.character(pi0))

	ret<-NULL
	ret$ssize<-ssize
	ret$power<-pwrMatrix
	ret$crit.vals<-critMatrix

	return(ret)		

}





