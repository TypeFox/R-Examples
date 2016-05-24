

ssize.Fvary<-function(X,beta,L=NULL,dn,a,b,fdr=0.05,power=0.8,pi0=0.95,maxN=20,cex.title=1.15,cex.legend=1){
	XTX<-t(X)%*%X
	B<-beta
	if (length(L)==0){L<-diag(length(B))}
	k<-length(L[1,])
	N<-maxN


	getAvgFcdf_varySigma<-function(c){
		intF<-function(rho){
			ncF<-pf(q=c,df1=k,df2=df,ncp=rho*lambda)*(1/(gamma(a)*(1/b)^a)*rho^(a-1)*exp(-rho*b))
			return(ncF)
		}
		val.int<-integrate(f=intF,lower=0,upper=Inf,abs.tol=1e-10)$value
		return(val.int)
	}


	FVarfun<-function(c){
		ratio<-fdr*(1-p)/((1-fdr)*p)
		dif<-abs((1-pf(q=c,df1=k,df2=df))/(1-getAvgFcdf_varySigma(c))-ratio)
		return(dif)	
	}

	crit<-NULL
	pwr2<-NULL
	ssize<-matrix(0,nrow=length(pi0),ncol=3)
	colnames(ssize)<-c("pi0","ssize","power")
	up.start<-100
	for(i in 1:length(pi0)){
		p<-pi0[i]; pwr.new<-0
        	up<-up.start
		for(n in 2:N){
			E<-t(L)%*%solve(t(X)%*%X)%*%L/n
			lambda<-t(t(L)%*%B)%*%solve(E)%*%(t(L)%*%B)
			df<-dn(n)
			ci<-optimize(f=FVarfun,interval=c(0,up))$min
			up<-ci    

			if((abs(ci-up.start)>=1)){pwr.new<-1-getAvgFcdf_varySigma(ci);crit.new<-ci}
			if((abs(ci-up.start)<1)&(pwr.new!=1)){pwr.new<-0;crit.new<-NA}
			crit<-c(crit,crit.new)	


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
	title(xlab="Sample size (n)",	ylab="Power")
	mtext(bquote("Average power vs. sample size with specified design matrix,"),
		cex=cex.title,padj=-2.35)
	mtext(bquote(paste("fdr=",.(round(fdr,4)),", and ",sigma[g]^2,"~IG(",.(round(a,4)),",",.(round(b,4)),")")),
		cex=cex.title,padj=-0.1)
	legend(x=N,y=0,xjust=1,yjust=0,col=1:i,pch=c(16,16,16),lty=1:length(pi0),
		legend=as.character(pi0),bg="white",title=expression(pi[0]),cex=cex.legend)

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
