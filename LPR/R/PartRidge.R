PartRidge<-function(x, y, lambda2=0, varset, standardize=TRUE, intercept=TRUE)
# compute the Ridge estimator: argmin (1/2n)|| Y-X\beta ||_2^2 + \lambda2/n || \beta ||_2^2
# lambda2: smoothing parameter. \lambda2=0 corresponds to Ordinary Least Squares
# varset: the set of variables not putting Ridge penalty
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

	penalty.factor<-rep(1,m)
	if(sum(varset)>0){
		penalty.factor[varset]<-0
	}
	xy<-t(x)%*%y/n
	xx<-t(x)%*%x/n
	diag(xx)<-diag(xx)+lambda2*penalty.factor
	#if(!is.singular.matrix(xx)){
	beta<-solve(xx,xy,tol=1e-64)
	#} else {
	#	LPR.fit<-glmnet(x=x, y=y, lambda=lambda2, alpha=0, standardize=FALSE, intercept=FALSE,
	#	penalty.factor=penalty.factor)
	#	beta<-LPR.fit$beta
	#}
	beta<-beta/normx  	
	object<-list()
	object$beta<-beta
	object$beta0<-mu-meanx%*%beta
	object$meanx<-meanx
	object$mu<-mu
	object$normx<-normx
	object$lambda2<-lambda2
	
	object    
}


  
  
  
