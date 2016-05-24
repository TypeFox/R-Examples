`penalized.pls.select` <-
function(X,y,M=NULL,ncomp,blocks){
    if (var(y)==0){
        coefficients=matrix(0,ncol(X),ncomp) # if y is constant, all coefficients are zero
        }
    else{

  X0<-X

  WW<-matrix(,ncol(X),ncomp)

  TT<-matrix(,nrow(X),ncomp)

  for (i in 1:ncomp){
    ww<-t(X)%*%y
    if (is.null(M)==FALSE){
        ww<-M%*%ww
    }
      # select optimal blocks
      score<-vector(length=max(blocks))
      for (k in 1:length(score)){
        score[k]=mean(ww[blocks==k]^2)
      }
      k.max<-which.max(score)
      w<-ww*0
      w[blocks==k.max]=ww[blocks==k.max]
      ww<-normalize.vector(w)
      #
      tt<-X%*%ww
      tt=normalize.vector(tt)
      WW[,i]=ww
      TT[,i]=tt
      X=X- tt%*%(t(tt)%*%X)
    }
    B=matrix(,ncol(X),ncomp)
        RR=matrix(t(TT)%*%X0%*%WW,nrow=ncomp) # the triangular matrix
        RR[row(RR)>col(RR)]<-0 # inserted for numerical stability
        LL=backsolve(RR,diag(ncomp))
        B=matrix(,ncol(X),ncomp)
        for (i in 1:ncomp) {
            Li=LL[1:i,1:i,drop=FALSE]
        B[,i]=WW[,1:i,drop=FALSE]%*%(Li%*%(t(TT[,1:i,drop=FALSE])%*%y))
    }
    coefficients=B
    }
  return(list(coefficients=coefficients))

}
