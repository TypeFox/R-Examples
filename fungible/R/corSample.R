##---------------------------------------------#
#  Generate sample correlation matrices
#  see Browne 1968
#  Browne, M. (1968). A comparison of factor analytic techniques. 
#  Psychometrika, 33(3):267-334.
#
# Output
#    cor.sample    a sampled correlation matrix
#    cov.sample    a sampled covariance matrix
#----------------------------------------------#

corSample<- function(R,n){
  Nvar<-ncol(R)
	Tmat<-matrix(0,Nvar,Nvar)
	Tmat[lower.tri(Tmat)]<-rnorm(n=(Nvar*(Nvar-1)/2))
	for(i in 1:Nvar){
# Note that Browne 1968 contains a typo for the df -- see the following for (n-i+1)
# Kshirsagar, A. (1959). Bartlett decomposition and wishart distribution. The Annals of Mathematical Statistics, 
#    30(1)239-241.
	 Tmat[i,i]<-sqrt(rchisq(n=1,df=(n-i+1)))
	}

	H<- Tmat %*% t(Tmat)
	
	Omega <-t(chol(R))
	A <-  Omega %*% H %*% t(Omega)
	S <- (1/n) * A
	Dmat<-diag(1/sqrt(diag(A)))
	R.samp<-Dmat%*%A%*%Dmat
	list(cor.sample=R.samp, cov.sample=S)
}


