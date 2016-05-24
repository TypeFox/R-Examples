`penalized.pls.kernel` <-
function(X,y,M=NULL,ncomp){
    if (var(y)==0){
        coefficients=matrix(0,ncol(X),ncomp) # if y is constant, the coefficients are zero
    }
    else{
  n=nrow(X)
  yhat=rep(0,n)
  if (is.null(M)==TRUE){
    K=X%*%t(X)
  }
  if (is.null(M)==FALSE){
    K=X%*%M%*%t(X)
  }
  UU=matrix(,n,ncomp)
  TT=matrix(,n,ncomp)
  for (i in 1:ncomp){
        uu=y-yhat
        uu=uu/sqrt(sum((K%*%uu)*uu))
        UU[,i]=uu
        if (i==1) {
            tt=K%*%uu
        }
        if (i>1){
            TTi=TT[,1:(i-1),drop=FALSE]
            tt = K%*%uu -  TTi %*% (t(TTi) %*% (K %*% uu))
            if (floor(i/5) == i/5) {
                tt = tt - TTi %*% (t(TTi) %*% tt)
            }
        }
        tt=normalize.vector(tt)
        TT[,i]=tt
        yhat=yhat + sum(tt*y) * tt
    }
    RR=t(TT)%*%K%*%UU
    RR[row(RR)>col(RR)]<-0 # inserted for numerical stability
    RR[row(RR)< col(RR) -1]<-0 # inserted for numerical stability
    LL=backsolve(RR,diag(ncomp))
    AA=matrix(,nrow(X),ncomp)
    for (i in 1:ncomp){
        Li=LL[1:i,1:i,drop=FALSE]
        AA[,i]= UU[,1:i,drop=FALSE]%*%(Li%*%(t(TT[,1:i,drop=FALSE])%*%y))
    }
    coefficients=t(X)%*%AA
    if (is.null(M)==FALSE){
    coefficients=M%*%coefficients   
}
}
  return(list(coefficients=coefficients))

}
