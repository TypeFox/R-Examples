S.LLR<-function (tt, h, Ker = Ker.norm,w=NULL,cv=FALSE)
{
 if (is.matrix(tt)) {
    if (ncol(tt)!=nrow(tt)) {
      if (ncol(tt)==1) {
         tt=as.vector(tt)
         tt=abs(outer(tt,tt, "-"))}
#      else stop("Error: incorrect arguments passed")
    }}
 else if (is.vector(tt))    tt=abs(outer(tt,tt, "-"))
 else stop("Error: incorrect arguments passed")
 if (cv)  diag(tt)=Inf
 k=Ker(tt/h)
 if (cv)   diag(k)=0
# S1=apply(k*tt,1,sum,na.rm=TRUE)
# S2=apply(k*(tt^2),1,sum,na.rm=TRUE)
 S1=rowSums(k*tt,na.rm=TRUE)
 S2=rowSums(k*(tt^2),na.rm=TRUE)
 b=k*(S2-tt*S1)
 if (cv)   diag(b)=0
if (is.null(w)) w<-rep(1,nrow(b))
 b<-sweep(b,1,w,FUN="*")
# res =b/apply(b,1,sum,na.rm=TRUE)
 res =b/rowSums(b,na.rm=TRUE)
return(res)}


