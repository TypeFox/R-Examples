design.matrix<-function(f1,f2){
	n<-length(f1)
	n1<-nlevels(as.factor(f1))
	n2<-nlevels(as.factor(f2))
	X0<-matrix(1,nrow=n,ncol=1)
	Xa<-matrix(0,nrow=n,ncol=n1)
	Xb<-matrix(0,nrow=n,ncol=n2)
	Xab<-matrix(0,nrow=n,ncol=n1*n2)
	for (i in 1:n1){
		ix<-which(f1==i)
		Xa[ix,i]<-1
	}
	for (i in 1:n2){
		ix<-which(f2==i)
		Xb[ix,i]<-1
	}
	list(X0=X0,Xa=Xa,Xb=Xb)
}
