# Loss aversion based c* calculation routines
# based on Esarey and Danneman, "A Quantitative Method for Substantive Robustness Assessment"
# Friday, December 18, 2009
# 
# Written by Justin Esarey, Rice University
#
#
# NOTE: These functions assume a t distributed posterior and are therefore appropriate
# for use after OLS regression models


cstar<-function(regobject, r){

	# translates regression results into a "loss point" using a simple kinked loss function
	# give it regression results, it returns a corresponding set of "loss points" t
	# if the point at which you begin making losses is less than t, accept evidence, else reject

	s<-sign(regobject$coefficients)
	b<-abs(regobject$coefficients)
	v<-sqrt(diag(vcov(regobject)))
	degfr<-regobject$df

	utility<-function(t,a,x,m,sigma,degfr){
		k<-as.numeric((x-t)>0)*2-1
		xx<-(x-m)/sigma;
		tt<-(t-m)/sigma;
		out<-a^(-k)*dt(xx, df=degfr)*(x-t);
		return(out)
	}

	maximand<-function(t,a,m,sigma,degfr){as.numeric(integrate(utility, m-8*sigma, m+8*sigma, rel.tol=.Machine$double.eps^0.5, t=t, a=a, m=m, sigma=sigma, degfr=degfr)[1])}
	
	out<-c()
	for(i in 1:length(b)){
		
		if(floor(i/10)==(i/10)){print(c("Iteration ",i))}
		if(b[i]!=0){out[i]<- sign(s[i])*uniroot(maximand, interval=c(-v[i]-b[i],v[i]+b[i]), tol=.Machine$double.eps^0.5, a=r, m=b[i], sigma=v[i], degfr=degfr)$root}
		else{out[i]<- 0}

	}

	for(i in 1:length(b)){

	if(sign(out[i])!=s[i]){out[i]<-0}

	}

	return(out)
	
}


cstari<-function(b,v,degfr,r){
	
	# "immediate" form of simple kink-slope statistic
	# give it a beta and a standard error and degrees of freedom and a gamma, it returns a "loss point" t
	
	s<-sign(b)
	b<-abs(b)	# reflect analysis into the positive domain (interpretation is the same)

	utility<-function(t,a,x,m,sigma,degfr){
		k<-as.numeric((x-t)>0)*2-1
		xx<-(x-m)/sigma;
		tt<-(t-m)/sigma;
		out<-a^(-k)*dt(xx, df=degfr)*(x-t);
		return(out)
	}

	maximand<-function(t,a,m,sigma,degfr){as.numeric(integrate(utility, b-8*sigma, b+8*sigma, rel.tol=.Machine$double.eps^0.5, t=t, a=a, m=m, sigma=sigma, degfr=degfr)[1])}
	
	out<-c()
	for(i in 1:length(b)){
		
		if(floor(i/10)==(i/10)){print(c("Iteration ",i))}
		if(b[i]!=0){out[i]<- sign(s[i])*uniroot(maximand, interval=c(-v[i]-b[i],v[i]+b[i]), tol=.Machine$double.eps^0.5, a=r, m=b[i], sigma=v[i], degfr=degfr)$root}
		else{out[i]<- 0}

	}

	for(i in 1:length(b)){

	if(sign(out[i])!=s[i]){out[i]<-0}

	}
	
	return(out)
	
}
