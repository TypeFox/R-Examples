.packageName <- "cheb"
`chebR` <-
function(a,b,tol=1e-15,relerr=0.0) {
m<-nrow(a); n<-ncol(a); ndim<-n+3; mdim<-m+1
if (n > m) stop("number of equations exceeds number of unknowns")
aa<-matrix(0,ndim,mdim); bb<-rep(0,mdim); xx<-rep(0,ndim)
aa[1:n,1:m]<-t(a); bb[1:m]<-b
rlist<-.Fortran("cheb",as.integer(m),as.integer(n),as.integer(m+1),as.integer(n+3),
	as.single(aa),bb=as.single(bb),as.single(tol),as.single(relerr),xx=as.single(xx),
	rank=as.integer(0),resmax=as.single(0.0),iter=as.integer(0),ocode=as.integer(0),PACKAGE="cheb")
return(list(coefs=rlist$xx[1:n],resids=rlist$bb[1:m],rank=rlist$rank,iter=rlist$iter,ocode=rlist$ocode))
}


