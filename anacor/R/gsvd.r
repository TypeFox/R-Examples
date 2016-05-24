gsvd <- function(r,p,q) {
	ps<-mfunc(p,function(x) ginvx(sqrt(x)))
	qs<-mfunc(q,function(x) ginvx(sqrt(x)))
	sv<-svd(ps%*%r%*%qs)
	return(list(gd=sv$d,gu=ps%*%sv$u,gv=qs%*%sv$v))
}