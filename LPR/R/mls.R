mls<-function(x,y,tau=0,standardize=TRUE,intercept=TRUE)
# compute the modified Least Squares estimator, tau=0 corresponds to the 
# Ordinary Least Squares estimator
{ 
	x<-as.matrix(x)
	n<-nrow(x)
	m<-ncol(x)
	one<-rep(1,n)
	if(intercept){
		meanx <- drop(one %*% x)/n
		x <- scale(x, meanx, FALSE)
		mu <- mean(y)
		y <- drop(y - mu)
	} else{
		meanx<-rep(0,m)
		mu<-0
		y<-drop(y)
	}
	if(standardize) {
		normx <- sqrt(drop(one %*% (x^2)))
		x <- scale(x, FALSE, normx)
		} else {
		normx <- rep(1, m)
	}

	if(abs(tau)<=1e-15){
		obj<-lm(y~.-1, data=data.frame(y=y, x=x))  
		beta<-coef(obj)
		beta[is.na(beta)]<-0		
	} else {
		s<-svd(x)
		Ux <- s$u
		Vx <- s$v
		Dx <- s$d/sqrt(n)
		indi<-Dx>tau
		D <-rep(0,length(Dx))
		D[indi]<-1/Dx[indi]
		beta<-1/sqrt(n)*Vx%*%diag(D,length(Dx),length(Dx))%*%t(Ux)%*%y		
	}

	beta<-scale(t(beta),FALSE,normx)
	object<-list()
	object$beta<-beta
	object$beta0<-mu-t(meanx)%*%t(beta) # intercept
	object$meanx<-meanx
	object$mu<-mu
	object$normx<-normx
	object$tau<-tau
	
	object  
}

