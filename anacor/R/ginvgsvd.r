ginvgsvd <- function(gs,p,ind) {
	x<-gs$gu; z<-gs$gv; gm<-gs$gd; gi<-gm[ind]
	iv1<-(ginvx(gm-gi)+ginvx(-(gm+gi)))/2
	iv2<-(ginvx(gm-gi)-ginvx(-(gm+gi)))/2
	n<-nrow(x); m<-nrow(z)
	a<-matrix(0,n+m,n+m)
	b<-tcrossprod(z%*%diag(iv2),x)
	a[1:n,n+(1:m)]<-t(b); a[n+(1:m),1:n]<-b
	a[n+(1:m),n+(1:m)]<-tcrossprod(z%*%diag(iv1),z)
	a[1:n,1:n]<-tcrossprod(x%*%diag(iv1),x)
	a[1:n,1:n]<-a[1:n,1:n]-(solve(p)-tcrossprod(x))/gi
	return(a)
}