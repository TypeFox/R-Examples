Pgpd<-function(y,x,N,Nexc,method,alpha,Z){
	Padth<-0.05
	z<-y[1:Nexc]
	t<-mean(y[Nexc:(Nexc+1)])
	z<-z-t
	z<-z^Z
	frac<-Nexc/N
	#Fitting the tail and computing the Pvalue
	if(method!="ML"&method!="MOM")stop("Wrong tail method")
	if(method=="ML"){
		parmhat<-c(1,0)
		myres<-optim(parmhat,nll,z=z)
		parmhat<-myres$par
		if(myres$convergence!=0)parmhat<-c(NA,NA)
		if(!any(is.na(parmhat))&parmhat[2]<0.5){
			a=parmhat[1];k=parmhat[2];n=length(z)
			mycov<-matrix(nrow=2,data=c(2*a^2/(1-k),a*(1-k),a*(1-k),(1-k)^2))
			mycov<-mycov/n
		}
		if(all(is.na(parmhat))){mycov<-matrix(nrow=2,ncol=2,data=NA)}
		if(!any(is.na(parmhat))&parmhat[2]>=0.5&parmhat[2]<1){
			mydt<-unlist(lapply(X=1:100,FUN=function(x) optim(c(1,0),nll,z=sample(z,replace=T))$par))
			mycov<-t(matrix(nrow=2,data=mydt))
		}
		if(any(is.na(parmhat))|parmhat[2]>=1){mycov<-matrix(nrow=2,ncol=2,data=NA)}
	}
	if(method=="MOM"){
		xbar<-mean(z)
		n<-length(z)
		s2<-var(z)
		a<-0.5*xbar*(xbar^2/s2 +1)
		k<-0.5*(xbar^2/s2 -1)
		parmhat<-c(a,k)
		mycov<-matrix(nrow=2,data=c(2*a^2*(1+6*k+12*k^2),a*(1+2*k)*(1+4*k+12*k^2),
		a*(1+2*k)*(1+4*k+12*k^2),(1+2*k)^2*(1+k+6*k^2)))
		mycov<-(mycov*(1+k)^2)/(n*(1+2*k)*(1+3*k)*(1+4*k))
	}
    	if(length(which(is.na(parmhat)))!=2){
        	a<-parmhat[1]
        	k<-parmhat[2]
		if(k==0){p<-exp(-z/a)} 
		if(k!=0){p<-(1-k*z/a)^(1/k)}
		if(k>0){p[z>a/k]<-0}
		if(!any(is.na(mycov))&nrow(mycov)==2){
			myeig<-eigen(mycov)
			V<-Re(myeig$vectors)
			D<-Re(diag(myeig$values))
			Q<-Re(matrix(ncol=2,data=rnorm(20000))%*%sqrt(abs(D))%*%t(V)+
			cbind(rep(parmhat[1],10000),rep(parmhat[2],10000)))
			k0vec<-which(abs(Q[,2])<.Machine$double.eps)
			kvec<-which(abs(Q[,2])>=.Machine$double.eps)
			mypci<-rep(1,10000)
			mypci[k0vec]<-exp(-x/Q[k0vec,1])
			mypci[kvec]<-(1-Q[kvec,2]*x/Q[kvec,1])^(1/Q[kvec,2])
			mypci[x>Q[,1]/Q[,2]&Q[,2]>0]<-0
			Phatci<-quantile(Re(mypci),c(alpha/2,1-alpha/2),na.rm=T)
		}
		if(nrow(mycov)==100){
			pos<-which(!is.na(mycov[,1]))
			if(length(pos)>50){
			mypci<-rep(1,100)
			k0vec<-which(abs(mycov[,2])<.Machine$double.eps)
                        kvec<-which(abs(mycov[,2])>=.Machine$double.eps)
			mypci[k0vec]<-exp(-x/mycov[k0vec,1])
			mypci[kvec]<-(1-mycov[kvec,2]*x/mycov[kvec,1])^(1/mycov[kvec,2])
			mypci[x>mycov[,1]/mycov[,2]&mycov[,2]>0]<-0
			Phatci<-quantile(Re(mypci),c(alpha/2,1-alpha/2),na.rm=T)
			}
		}
		if(all(is.na(mycov)))Phatci<-c(NA,NA)
		p<-t(sort(p))
		n<-length(p)
		i<-1:n
		#Cramer-von Mises statistic
		#W2<-sum((p-((2*i-1)/(2*n)))^2)+1/(12*n)
		#Anderson Darling statistic
		A2<- -n-(1/n)*((2*i-1)%*%t(log(p)+log(1-p[n+1-i])))
		#load tables
		ktable<-matrix(nrow=10,ncol=1,data=c(0.9,0.5,0.2,0.1,0.0,-0.1,
			-0.2,-0.3,-0.4,-0.5))
		ptable<-matrix(nrow=1,ncol=8,data=c(0.5,0.25,0.1,0.05,0.025,0.01,
			0.005,0.001))
		A2table<-t(matrix(nrow=8,ncol=10,data=c(0.339,0.471,0.641,0.771,0.905,1.086,
			1.226,1.559,
			0.356,0.499,0.685,0.830,0.978,1.180,1.336,1.707,
			0.376,0.534,0.741,0.903,1.069,1.296,1.471,1.893,
			0.386,0.550,0.766,0.935,1.110,1.348,1.532,1.966,
			0.397,0.569,0.796,0.974,1.158,1.409,1.603,2.064,
			0.410,0.591,0.831,1.020,1.215,1.481,1.687,2.176,			
			0.426,0.617,0.873,1.074,1.283,1.567,1.788,2.314,
			0.445,0.649,0.924,1.140,1.365,1.672,1.909,2.475,			
			0.468,0.688,0.985,1.221,1.465,1.799,2.058,2.674,
			0.496,0.735,1.061,1.321,1.590,1.958,2.243,2.922)))
		myk<-max(0.5,k)
		myp<-apply(A2table,2,FUN=function(x) interp1(ktable[,1],x,myk,
			"linear",extrap=TRUE))
		Pad<-max(min(interp1(myp,t(ptable),A2,"linear",extrap=TRUE),1),0)
        	if(Pad>Padth){
			x<-(x-t)^Z
			if(abs(k)<.Machine$double.eps)Phat<-exp(-x/a)  
			if(abs(k)>=.Machine$double.eps)Phat<-(1-k*x/a)^(1/k)   
			#if(x>(a/k)&k>0)Phat<-mean(mypci,na.rm=T)
			if(x>(a/k)&k>0)Phat<-0
			Phat = frac*Phat
			Phatci = frac*Phatci;
		}else{Phat<-NA;Phatci<-c(NA,NA);k<-NA}
	}else{Phat<-NA;Phatci<-c(NA,NA);k<-NA}
	return(list(Phat,Phatci,k))
}
