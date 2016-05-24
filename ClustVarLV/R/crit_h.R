crit_h <-
function(method,X12,EXTr,Xr,EXTu,Xu12)
### clustering criterion / hierarchy
{
  n<-nrow(X12)
  pk<-ncol(X12)  
  
  if (method == 1) {
    if((EXTr==0)&(EXTu==0)) {
      xnew<- X12  
      if (dim(xnew)[1]>dim(xnew)[2]) vp = eigen( t(xnew) %*% xnew)
      else vp = eigen( xnew %*% t(xnew))
      crit = vp$values[1] /(n-1)  
    }   
    if((EXTr==1)&(EXTu==0)) {
      xnew<- t(Xr)%*%X12  
      if (dim(xnew)[1]>dim(xnew)[2]) vp = eigen( t(xnew) %*% xnew)
      else vp = eigen( xnew %*% t(xnew))
      crit = vp$values[1] /(n-1)  
    }
    if((EXTr==0)&(EXTu==1))  {
      P<-X12 %*% Xu12
      B<-t(X12)%*%P
      vp <- eigen(t(B) %*% B)
      alpha2<-eigen(t(P)%*%P)$values[1]
      crit= vp$values[1]/((n-1)*alpha2)
    }
  }
  
  if (method==2) {
    if ((EXTu==0)&(EXTr==0)){          
      xbar = X12 %*% matrix(1,pk,1) /pk
      #crit = pk * 1/(n-1) * t(xbar) %*% xbar         # version RSA      ck = xbark
      crit = pk * sqrt(1/(n-1) * t(xbar) %*% xbar)    # version CommStat ck normalized
    }
    if ((EXTu==0)&(EXTr==1)){                              
      xbar = X12 %*% matrix(1,pk,1) /pk                
      pxbar=sqrt (t(xbar)%*%Xr%*%t(Xr)%*%xbar)
      crit= pk * 1/(n-1) *pxbar 
    }
    if ((EXTu==1)&(EXTr==0)){
      P=X12 %*% Xu12
      alpha2<-sum(diag(t(P)%*%P))
      xbar = X12 %*% matrix(1,pk,1) /pk       
      pxbar=sqrt (t(xbar)%*%P%*%t(P)%*%xbar)
                 if (is.nan(pxbar)) { print(Xu12)}
      crit<- pxbar*pk/sqrt(n-1)
      crit<-crit/sqrt(alpha2)
    }
  }
 
return(crit)
}
