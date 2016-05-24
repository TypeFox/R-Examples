"norm.fdata"<-function(fdataobj,metric=metric.lp,...){
if (!inherits(fdataobj,"fdata")) stop("No fdata class")
if (is.vector(fdataobj$data))    fdataobj$data=matrix(fdataobj$data,nrow=1)
z0<-matrix(0,ncol=ncol(fdataobj),nrow=1)
z0<-fdata(z0,fdataobj[["argvals"]],fdataobj[["rangeval"]],fdataobj[["names"]])
n.lp<-metric(fdataobj,z0,...)
n.lp
}


"norm.fd"<-function(fdobj){
if (is.fd(fdobj)) rng<- fdobj[[2]]$rangeval
else if (is.basis(fdobj)) rng<- fdobj$rangeval
else stop("Non fd or basis class")
sqrt(inprod(fdobj,fdobj))#/diff(rng))
}


