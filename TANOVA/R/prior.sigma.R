prior.sigma<-function(Q,f1,f2,tp=0){
	n<-dim(Q)[1]
	n1<-nlevels(as.factor(f1))
	if (length(tp)>1){
		n3<-nlevels(as.factor(tp))
		if (length(f2)>1){
			n2<-nlevels(as.factor(f2))
			df<-length(f1)-n1*n2*n3
		}
		else{
			df<-length(f1)-n1*n3
		}
	}
	if (length(tp)==1){
		if (length(f2)>1){
			n2<-nlevels(as.factor(f2))
			df<-length(f1)-n1*n2
		}
		else{
			df<-length(f1)-n1
		}
	}
	sg<-apply(Q,1,sum)/df
	zg<-log(sg)
	eg<-zg-digamma(df/2)+log(df/2)
	e<-mean(eg)
	x<-mean((eg-e)^2*n/(n-1)-trigamma(df/2))
	if (x<=0){
		v0<-10^5
	}
	else{
		v0<-2*trigammaInverse(x)
	}
	s0<-exp(e+digamma(v0/2)-log(v0/2))
return(list(s0=s0,v0=v0,df=df))
}