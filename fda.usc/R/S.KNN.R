S.KNN<-function(tt,h=NULL,Ker=Ker.unif,w=NULL,cv=FALSE)      {
if (is.matrix(tt)) {
    if (ncol(tt)!=nrow(tt)) {
      if (ncol(tt)==1) {
         tt=as.vector(tt)
         tt=abs(outer(tt,tt, "-"))}
#      else stop("Error: incorrect arguments passed")
    }}
 else if (is.vector(tt))    tt=abs(outer(tt,tt, "-"))
 else stop("Error: incorrect arguments passed")
numgr=ncol(tt)
if (is.null(h)) h=floor(quantile(1:numgr,probs=0.05,na.rm=TRUE,type=4))
else if (h<=0 ) stop("Error: incorrect knn value")
tol=1e-19
tol=diff(range(tt)*tol)
tol=1e-19
if (cv) diag(tt)=Inf
vec=apply(tt,1,quantile,probs=((h)/numgr),type=4)+tol
rr=sweep(tt,1,vec,"/")
rr=Ker(rr)
#if (cv) diag(rr)=0
if (!is.null(w)){ #w<-rep(1,ncol(rr))
  rr<-sweep(rr,2,w,FUN="*") }  ## antes un 2
#print(colSums(rr,na.rm=TRUE))
rr=rr/rowSums(rr,na.rm=TRUE)
return(rr)
}
