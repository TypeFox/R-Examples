mat_init <-
function(X,EXTr,Xr,EXTu,Xu,K) 
{
  
  if ((EXTr==0)|(EXTu==0)) { 
    n<-nrow(X)
    comp<-matrix(0,nrow=n,ncol=K)
    rownames(comp)<-rownames(X)
    colnames(comp)<-paste("Comp",c(1:K),sep="")
    a<-NULL
    u<-NULL
  }
  if (EXTr==1) {
    q<-ncol(Xr)
    a <-matrix(0,nrow=q,ncol=K)
    rownames(a)<-colnames(Xr)
    colnames(a)<-paste("Comp",c(1:K),sep="")
    u<-NULL
  }
  if (EXTu==1) {
    m<-ncol(Xu)
    #u<-matrix(0,nrow=m,ncol=K)
    u<-matrix(runif(m*K),nrow=m,ncol=K)
    for (k in 1:K) {u[,k]<-u[,k]/as.numeric(sqrt(t(u[,k])%*%u[,k]))}
    rownames(u)<-colnames(Xu)
    colnames(u)<-paste("Comp",c(1:K),sep="")
    a<-NULL
  }     
mat<-list(comp=comp,a=a,u=u)  
return(mat)
}
