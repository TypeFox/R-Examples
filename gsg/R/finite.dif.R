## finite differences first and secon partial derivatives
## because the relationship of population mean fitness
## to population mean fitness will be quite smooth relative
## the the individual fitness function, this should be
## perfectly precise, requiring few function evaluations
## get the partial derivatives, and without necessarily
## demanding extremely fine precision of the multivariate 
## integral of of fitness
finite.dif<-function(ff,x,h,quiet=TRUE,...){
	count<-0
	if(quiet==FALSE){
		print(paste(count,"of",2*length(x)+(length(x)^2-length(x))/2+1,
		     "function evaluations completed."))
	}
	grad.est<-array(dim=length(x))
	hessian.est<-matrix(NA,length(x),length(x))
	ffx<-ff(x,...)
	count<-count+1
	if(quiet==FALSE){
		print(paste(count,"of",2*length(x)+(length(x)^2-length(x))/2+1,
		     "function evaluations completed."))
	}
	ffxx<-matrix(NA,length(x),length(x))
	for(i in 1:length(x)){
		x1<-x2<-x
		x1[i]<-x[i]+h[i]
		x2[i]<-x[i]-h[i]
		ffx1<-ff(x1,...)
	count<-count+1
	if(quiet==FALSE){
		print(paste(count,"of",2*length(x)+(length(x)^2-length(x))/2+1,
		     "function evaluations completed."))
	}
		ffxx[i,i]<-ffx1
		ffx2<-ff(x2,...)
	count<-count+1
	if(quiet==FALSE){
		print(paste(count,"of",2*length(x)+(length(x)^2-length(x))/2+1,
		      "function evaluations completed."))
	}
		grad.est[i]<-(ffx1-ffx2)/(2*h[i])
		hessian.est[i,i]<-((ffx1-ffx) - (ffx-ffx2))/h[i]^2
	}
	for(i in 1:(length(x)-1)){
		for(j in (i+1):length(x)){
			x1<-x
			x1[i]<-x[i]+h[i]
			x1[j]<-x[j]+h[j]
			ffxx[i,j]<-ff(x1,...)
	count<-count+1
	if(quiet==FALSE){
		print(paste(count,"of",2*length(x)+(length(x)^2-length(x))/2+1,
		     "function evaluations completed."))
	}
		}
	}
	for(i in 1:(length(x)-1)){
		for(j in (i+1):length(x)){
			hessian.est[i,j]<-hessian.est[j,i]<-
			    ((ffxx[i,j]-ffxx[j,j])/h[i] - (ffxx[i,i]-ffx)/h[i])/h[j]
		}
	}
	list(grad.est=grad.est,hessian.est=hessian.est,mu=ffx)
}
