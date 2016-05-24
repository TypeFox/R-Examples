
computeAuxVariables_update <- function(model) {
	
	aux <- covMatrix(X=model@X, object=model@covariance, noise.var=model@noise.var)
	C <- aux[[1]]
	
	T <- chol(C)
	t.T<-t(T)
	
	M <- backsolve(t.T, model@F, upper.tri = FALSE)
	z <- backsolve(t.T, model@y-model@F%*%as.matrix(model@trend.coef), upper.tri=FALSE) 	
	
	model@T <- T
	model@z <- as.numeric(z)
	model@M <- M
	return(model)
} 


#computeAuxVariables.update2 <- function(model,oldchol=NULL,newX=NULL,X=NULL) {
#	
#	n<-model@n
#	
#	if(n>20){
#		if(is.null(oldchol)) oldchol<-t(model@T)
#		if(is.null(X)) X <- as.matrix(model@X)
#		if(is.null(newX)) newX <- matrix(X[model@n,],ncol=model@d)
#	
#		newcov=covMat1Mat2(X1=newX, X2=X, object=model@covariance,nugget.flag=model@covariance@nugget.flag)
#		newcov<-as.numeric(newcov)	
#		t.T <- chol.update(newcov=newcov,old.chol=oldchol,return.all=TRUE)
#	}else{
#		aux <- covMatrix(X=model@X, object=model@covariance, noise.var=model@noise.var)
#		C <- aux[[1]]		
#		T <- chol(C)
#		t.T<-t(T)
#	}
#	
#	M <- backsolve.update(t.T, model@F,model@M)
#	z <- backsolve.update(t.T, model@y-model@F%*%model@trend.coef, matrix(model@z,ncol=1)) 	
#	
#	model@T <- t(t.T)
#	model@z <- as.numeric(z)
#	model@M <- M
#	return(model)
#} 
#
#backsolve.update<-function(U,V,prev.sol){
#	
#	#solve U x = V
#	#U is a lower triangular matrix
#	#prev.sol is solution of the backsolve at step n-1
#	#sum(U[n,i]) prevsol[i,] + U[n,n] . x = V[n]
#		
#	n<-nrow(U)
#	u1<-U[n,]
#	v1<-V[n,]
#	
#	u1.trunc <- u1[-n]
#	u1.n <- u1[n]
#	tpvect<-crossprod(u1.trunc,prev.sol)
#	
#	sol.n <- as.numeric((v1-tpvect)/u1.n)
#	return(rbind(prev.sol,sol.n))
#	
#}
#
#chol.update<-function(mat=NULL,newcov=NULL,old.chol,return.all=TRUE){
#	
#	#old.chol is a n x n lower triangular matrix
#	#you can provide either mat or newcov
#	#mat is a n+1 x n+1 matrix
#	#newcov is a n+1 vector
#	
#	if(is.null(newcov)){
#		if(is.null(mat)) return("error")
#		n<-ncol(mat)
#		newcov<-mat[,n]
#	}
#	p<-length(newcov)
#	
#	b<-newcov[-p]
#	w<-newcov[p]
#	
#	x<-backsolve(old.chol,b,upper.tri=FALSE)
#	#x<-rep(0,times=(p-1))
#	
#	z<-as.numeric(sqrt(w-sum(x*x)))
#	
#	if(!return.all) return(c(x,z))
#	newmat<-cbind(old.chol,rep(0,times=p-1))
#	newmat<-rbind(newmat,c(x,z))
#	return(newmat)
#	
#}
