#proj <- function(x, ...) UseMethod("proj") 

proj.dir<-function(data,f1,f2,time.course,type,eb=FALSE,df=0,...){
	n<-dim(data)[1]/time.course
	a<-matrix(nrow=n,ncol=time.course)
	if (eb==TRUE){
		pr<-prior.SIGMA(data,f1,f2,time.course)
	}
	if (length(f2)==1){
		e<-matrix(1,nrow=length(f1),ncol=1)
		for (i in 1:n){
			ix<-c(((i-1)*time.course+1):(i*time.course))
			x<-data[ix,]
			if (eb==TRUE){
				sigma<-(sigma.hat(x,f1,f2)*pr$df+pr$LAMBDA*pr$v0)/(pr$df+pr$v0)
			}
			else{
				sigma<-sigma.hat(x,f1,f2)
			}
			a[i,]<-Re(eigen((x%*%t(x)-1/length(f1)*x%*%e%*%t(e)%*%t(x))%*%ginv(sigma))$vector[,1])
			a[i,]<-sign(a[i,1])*a[i,]
		}
		return(a)
	}
	
	if (length(f2)>1){
		if (type==1){
			m<-ls.estimate(data,f1,f2)
			z<-data-m$Mab
		}
		if (type==2){
			m<-ls.estimate(data,f1,f2)
			z<-data-m$M0
		}
		if (type==3){
			m<-ls.estimate(data,f1,f2)
			z<-data-m$Mb
		}
		if (type==4){
			m<-ls.estimate(data,f2,f1)
			z<-data-m$Mb
		}
		for (i in 1:n){
			ix<-c(((i-1)*time.course+1):(i*time.course))
			x<-z[ix,]
			if (eb==TRUE){
				sigma<-(sigma.hat(x,f1,f2)*pr$df+pr$LAMBDA*pr$v0)/(pr$df+pr$v0)
			}
			else{
				sigma<-sigma.hat(x,f1,f2)
			}
			if (df==0){
				a[i,]<-Re(eigen(x%*%t(x)%*%ginv(sigma))$vector[,1])
				a[i,]<-sign(a[i,1])*a[i,]
			}
			else {
				N<-ns(c(1:time.course),df=min(df,time.course),intercept=TRUE)
				theta<-Re(eigen(ginv(t(N)%*%(N))%*%t(N)%*%x%*%t(x)%*%ginv(sigma)%*%N)$vector[,1])
				a[i,]<-N%*%theta
				a[i,]<-sign(a[i,1])*a[i,]
a[i,]<-a[i,]/sqrt(a[i,]%*%a[i,])
}
}
return(a)
}
}
