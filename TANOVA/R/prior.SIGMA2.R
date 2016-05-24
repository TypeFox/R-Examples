prior.SIGMA<-function(data,f1,f2,time.course){
	n<-dim(data)[1]/time.course
	n1<-nlevels(as.factor(f1))
	if (length(f2)>1){
		n2<-nlevels(as.factor(f2))
		dg<-length(f1)-n1*n2
	}
	if (length(f2)==1){
		dg<-length(f1)-n1
	}
	s<-matrix(nrow=time.course,ncol=n)
	z<-matrix(0,nrow=time.course,ncol=time.course)
	for (i in 1:n){
		ix<-c(((i-1)*time.course+1):(i*time.course))
		temp<-sigma.hat(data[ix,],f1,f2)
		s[,i]<-diag(temp)
		z<-z+temp
	}
	v<-vector(length=time.course)
	for (i in 1:time.course){
		sg<-s[i,]
		zg<-log(sg)
		eg<-zg-digamma(dg/2)+log(dg/2)
		e<-mean(eg)
		x<-mean((eg-e)^2*n/(n-1)-trigamma(dg/2))
		if (x<=0){
			v[i]<-10^5
		}
		else{
			v[i]<-2*trigammaInverse(x)
		}
	}
	v0<-max(mean(v),time.course+6)
	LAMBDA<-1/v0*(v0-time.course-1)*z/(n-1)
	v0<-mean(v)
	return(list(LAMBDA=LAMBDA,v0=v0,df=dg))
}