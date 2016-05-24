###############################################################################
# TANOVA: test.R
# 
# TODO: Add comment
#
# Author: Weihong
# Mar 20, 2009, 2009
###############################################################################
ls.estimate<-function(data,f1,f2){
	n<-dim(data)[1]
	p<-length(f1)
	m<-design.matrix(f1,f2)
	Mab<-matrix(nrow=n,ncol=p)
	Mb<-matrix(nrow=n,ncol=p)
	M0<-matrix(nrow=n,ncol=p)
	X<-cbind(m$X0,m$Xa,m$Xb)
	temp<-X%*%ginv(t(X)%*%(X))%*%t(X)
	Mab<-t(temp%*%t(data))
	X<-cbind(m$X0,m$Xb)
	temp<-X%*%ginv(t(X)%*%(X))%*%t(X)
	Mb<-t(temp%*%t(data))
	X<-m$X0
	temp<-X%*%ginv(t(X)%*%(X))%*%t(X)
	M0<-t(temp%*%t(data))
	return (list(Mab=Mab,Mb=Mb,M0=M0))
}
