consol_affect <-
function(method,X,Xr,Xu,EXTr,EXTu,comp,a,u)
{   
  p<-ncol(X)
  gtmp<-rep(0,p)
  n<-nrow(X)
 
  if (method == 1){
    if ((EXTu==0)&(EXTr==0)) {
      comp<-comp/matrix(sqrt(diag(t(comp)%*%comp)),nrow=n,ncol=ncol(comp),byrow=T)  #normalization de comp
      cova = t(comp) %*% X/(n-1)
      covaC<-  cova^2
      for (j in (1:p)) {   
        vec=covaC[,j] 
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    } 
    if ((EXTu==0)&(EXTr==1)) {
      cova = t(comp) %*% X/(n-1)
      covaC <- cova^2
      for (j in (1:p)) {
        vec=covaC[,j]
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    }  
    if ((EXTu==1)&(EXTr==0)){          
      cova = t(comp) %*% X/(n-1)
      covaC <- cova^2
      for (j in (1:p)) {
        vec=covaC[,j]
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    }   
  }
  
  if (method==2){
    if ((EXTu==0)&(EXTr==0)) {
      comp<-comp/matrix(sqrt(diag(t(comp)%*%comp)),nrow=n,ncol=ncol(comp),byrow=T)  #normalization de comp
      cova = t(comp) %*% X/(n-1)
      covaC<-  cova
      for (j in (1:p)) {
        vec=cova[,j]
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    } 
    if ((EXTu==0)&(EXTr==1)) {
      cova = t(comp) %*% X/(n-1)
      covaC <- cova
      for (j in (1:p)) {
        vec=cova[,j]
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    }  
    if ((EXTu==1)&(EXTr==0)){          
      for (j in 1:p) {
        Pj<- as.matrix(X[,j])%*%Xu[j,]
        critind=diag(t(u)%*%t(Pj)%*%comp)  
        vec=critind
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]                    
      }
    }   
  } 
  
  return(gtmp)
}
