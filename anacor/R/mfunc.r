mfunc<-function(a,fn=sqrt) {
	e<-eigen(a)
  y<-e$vectors
  v<-e$values
	return(tcrossprod(y%*%diag(fn(v)),y))
}