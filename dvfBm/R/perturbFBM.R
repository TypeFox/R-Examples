perturbFBM<-function(n,H,C=1,type="no",SNR=NULL,plot=FALSE){

################################################################################
## provides a contamindated sample path of a fBm with Hurst parameter 
## H and scaling coefficient C, at times i=1,...,n
## see Achard and Coeurjolly (2009) for a detail of the contamination parameters
################################################################################

	z<-circFBM(n,H,plot=FALSE)*n^H
	switch(type,
		"no"={
			zz<-z
		},
		"B0"={ 
			s2B<-10^(-SNR/10)*C^2
			zz<-z+cumsum(rnorm(n,0,sqrt(s2B)))
		},
		"B1"={	
			s2B<-.5*10^(-SNR/10)*C^2
			zz<-z+rnorm(n,0,sqrt(s2B))
		},
		"AO"={
			p<-.005;ind<-which(rbinom(n-1,1,p)==1)
			if (length(ind)==0) ind<-sample(1:(n-1),1)
			nb<-length(ind)
			tmp<-diff(z)
			s2Pert<- 10^(-SNR/10)*C^2
			pert<-sqrt(s2Pert)*rnorm(nb)
			tmp[ind]<-tmp[ind]+pert
			zz<-c(0,cumsum(tmp))			
           	}
		)
	if (plot) {
		txtMain<-paste("n=",n,", H=",H,sep="")
		if (type=="no") {
			par(mfrow=c(2,1))
			plot(z,type="l",main=txtMain,xlab="Time",ylab="FBM");plot(diff(z),type="l",main=txtMain,xlab="Time",ylab="FGN")
		}
		if (type %in% c("B0","B1")) {
			txtMain2<-paste("n=",n,", H=",H,", SNR=",SNR,sep="")
			par(mfrow=c(2,2))
			plot(z,type="l",ylim=range(c(z,zz)),main=txtMain,xlab="Time",ylab="FBM");
			plot(zz,type="l",ylim=range(c(z,zz)),col='red',main=txtMain2,xlab="Time",ylab="noisy FBM")
			plot(diff(z),type="l",ylim=range(c(diff(z),diff(zz))),main=txtMain,xlab="Time",ylab="FGN")
			plot(diff(zz),type="l",col='red',ylim=range(c(diff(z),diff(zz))),main=txtMain2,xlab="Time",ylab="noisy FGN")
		}
		if (type %in% c("AO")) {
			txtMain2<-paste("n=",n,", H=",H,", SNR=",SNR,sep="")
			par(mfrow=c(2,1))
			plot(z,type="l",ylim=range(c(z,zz)),main=txtMain2,xlab="Time",ylab="FBM + outliers");
			segments(ind,z[ind],ind,zz[ind])
			points(ind,zz[ind],type="p")
			plot(diff(z),type="l",ylim=range(c(diff(z),diff(zz))),main=txtMain2,xlab="Time",ylab="FGN + outliers")
			segments(ind,diff(z)[ind],ind,diff(zz)[ind])
			points(ind,diff(zz)[ind],type="p")
		}
		par(mfrow=c(1,1))
	}
	zz
}
