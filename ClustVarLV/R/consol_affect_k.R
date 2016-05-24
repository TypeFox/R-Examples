consol_affect_k <-
function(method,X,Xr,Xu,EXTr,EXTu,comp,a,u,rlevel)
{   
  p<-ncol(X)
  gtmp<-rep(0,p)
  n<-nrow(X)
  comp<-comp/matrix(sqrt(diag(t(comp)%*%comp)),nrow=n,ncol=ncol(comp),byrow=T)  #normalization de comp
  if (method == 1){
    if ((EXTu==0)&(EXTr==0)) {

      cova = t(comp) %*% X/(n-1)
      covaC<-  cova^2
      for (j in (1:p)) {
        
        vec=c(covaC[,j],rlevel^2*var(X[,j])/(n-1))
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    } 
    if ((EXTu==0)&(EXTr==1)) {
      cova = t(comp) %*% X/(n-1)
      covaC <- cova^2
      for (j in (1:p)) {
     
        vec=c(covaC[,j],rlevel^2*var(X[,j])/(n-1))
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    }  
    if ((EXTu==1)&(EXTr==0)){          
      cova = t(comp) %*% X/(n-1)
      covaC <- cova^2
      for (j in (1:p)) {
      
        vec=c(covaC[,j],rlevel^2*var(X[,j])/(n-1))
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    }   
  }
  
  if (method==2){
    if ((EXTu==0)&(EXTr==0)) {
    
      cova = t(comp) %*% X/(n-1)
      covaC<-  cova
      for (j in (1:p)) {
      
        vec=c(cova[,j],rlevel*sqrt(var(X[,j])/(n-1)))
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    } 
    if ((EXTu==0)&(EXTr==1)) {
     
      cova = t(comp) %*% X/(n-1)
      covaC <- cova
      for (j in (1:p)) {
        vec=c(cova[,j],rlevel*sqrt(var(X[,j])/(n-1)))
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
    
      }
    }  
    if ((EXTu==1)&(EXTr==0)){          
      for (j in 1:p) {
        
        Pj<- as.matrix(X[,j])%*%Xu[j,]
        critind=diag(t(u)%*%t(Pj)%*%comp)  
  
        vec=c(critind,rlevel*sqrt(var(X[,j])/(n-1)))
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]                    
      }
    }   
  } 
  
  return(gtmp)
}
