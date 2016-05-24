crit_init <-
function(method,X,EXTr,Xr,EXTu,Xu)
### initial clustering criterion (each variable = one group) / hierarchy
{
  n<-nrow(X)
  p<-ncol(X)
 
  # method 1
  if (method==1) {
    if ((EXTu==0)&(EXTr==0)){ crit<- apply(X^2/(n-1), 2, sum) } 
    if ((EXTu==0)&(EXTr==1)){ 
      XrX<- t(Xr)%*%X
      crit<- apply(XrX^2/(n-1), 2, sum)
    }
    if ((EXTu==1)&(EXTr==0)){
      crit=c()
      for (i in 1:p) {
        critk<-var(X[,i])
        crit=c(crit,critk)
      }
    }
  }
  # method 2                 
  if (method==2) {
     if ((EXTu==0)&(EXTr==0)){ 
      #crit<-apply(X,2,var)    # version RSA        ck=xbark
      crit<-apply(X,2,sd)      # version CommStat   ck normalized
    } 
    if ((EXTu==0)&(EXTr==1)){ 
      px<-sqrt (diag(t(X)%*%Xr%*%t(Xr)%*%X))
      crit<- px/(n-1)
    } 
    if ((EXTu==1)&(EXTr==0)){
      crit=c()
      for (i in 1:p) {
        critk<- sqrt((t(X[,i])%*%X[,i])/(n-1))
        crit=c(crit,critk)
      }
    }
  }
  
  return(crit)
}
