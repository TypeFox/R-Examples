S.basis=function(tt,basis,lambda=0,Lfdobj=vec2Lfd(c(0,0)),w=NULL,...){
    phi=getbasismatrix(tt,basis,...)
    np<-length(tt)
    if (is.null(w)) w<-rep(1,np)
    if (!is.matrix(w)) w<-diag(w)
    if (lambda!=0) {
       R=getbasispenalty(basis,Lfdobj)
       S=phi%*%solve(t(phi)%*%w%*%phi+lambda*R)%*%t(phi)%*%w}
   else {S=phi%*%solve(t(phi)%*%w%*%phi)%*%t(phi)%*%w}
    return(S)
}
