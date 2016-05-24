S.NW<-function (tt, h=NULL, Ker = Ker.norm,w=NULL,cv=FALSE) {
 if (is.matrix(tt)) {
    if (ncol(tt)!=nrow(tt)) {
      if (ncol(tt)==1) {
         tt=as.vector(tt)
         tt=abs(outer(tt,tt, "-"))}
      #else stop("Error: incorrect arguments passed")
    }}
 else if (is.vector(tt))    tt=abs(outer(tt,tt, "-"))
 else stop("Error: incorrect arguments passed")
 if (is.null(h)) {
    h=quantile(tt,probs=0.15,na.rm=TRUE)
    cat("h=");print(h)
    }
  if (cv)  diag(tt)=Inf
  tt2<-as.matrix(sweep(tt,1,h,FUN="/"))
  k<-Ker(tt2)
  if (is.null(w)) w<-rep(1,len=ncol(tt))
  k1<-sweep(k,2,w,FUN="*")
#  S =k1/apply(k1,1,sum)
  rw<-rowSums(k1)
  rw[rw==0]<-1
  S =k1/rw
return(S)
}
